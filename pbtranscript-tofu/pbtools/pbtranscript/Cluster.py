"""Define class `Cluster` and `ClusterException`."""
#!/usr/bin/env python
import os
import os.path as op
import logging
import cPickle
from pbcore.io import FastaReader
from pbtools.pbtranscript.PBTranscriptException import PBTranscriptException
from pbtools.pbtranscript.io.FastaSplitter import splitFasta
from pbtools.pbtranscript.Utils import realpath, ln
from pbtools.pbtranscript.Polish import Polish
#from pbtools.pbtranscript.ice.ProbModel import ProbFromModel, ProbFromQV
from pbtools.pbtranscript.c_Prob import ProbFromModel, ProbFromQV # this replaces the above for speed optimization
from pbtools.pbtranscript.io.Summary import ClusterSummary
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ice.IceInit import IceInit
from pbtools.pbtranscript.ice.IceIterative import IceIterative
from pbtools.pbtranscript.ice.IceUtils import clean_up_after_ICE
from pbtools.pbtranscript.__init__ import get_version


class ClusterException(PBTranscriptException):
    """
    Exception class for Classifier.
    """
    def __init__(self, msg):
        PBTranscriptException.__init__(self, "cluster", msg)


class Cluster(IceFiles):
    """
    An object of `Cluster` calls the ICE algorithm to
    generate consensus isoforms.
    """
    def __init__(self, root_dir, flnc_fa, nfl_fa,
                 bas_fofn, ccs_fofn, out_fa,
                 sge_opts, ice_opts,
                 hq_isoforms_fa=None, hq_isoforms_fq=None,
                 lq_isoforms_fa=None, lq_isoforms_fq=None,
                 report_fn=None, summary_fn=None):
        super(Cluster, self).__init__(prog_name="Cluster",
                root_dir=root_dir, bas_fofn=bas_fofn, ccs_fofn=ccs_fofn)

        self.flnc_fa, self.nfl_fa, self.ccs_fofn = self._validateInputs(
            flnc_fa, nfl_fa, ccs_fofn)

        self.root_dir, self.out_fa = self._validateOutputs(
            root_dir, out_fa)

        self.hq_isoforms_fa = hq_isoforms_fa
        self.hq_isoforms_fq = hq_isoforms_fq
        self.lq_isoforms_fa = lq_isoforms_fa
        self.lq_isoforms_fq = lq_isoforms_fq

        self.sge_opts = sge_opts  # SGE, CPU options and etc
        self.ice_opts = ice_opts  # The ICE clutering algorithm options

        self.sanity_check()

        self._probqv = None     # probability & quality value

        self._flnc_splitted_fas = []  # split flnc_fa into smaller files.
        self._nflncSplittedFas = []  # split nfl_fa into smaller files.
        self._logConfigs()      # Log configurations

        self.iceinit = None
        self.icec = None
        self.iceq = None
        self.pol = None

        self.report_fn = realpath(report_fn) if report_fn is not None \
                else op.join(self.root_dir, "cluster_report.csv")
        self.summary_fn = realpath(summary_fn) if summary_fn is not None \
                else op.join(self.root_dir, "cluster_summary.txt")

        self.summary = ClusterSummary()

        self.add_log("Finishing create Cluster Object.", level=logging.INFO)

    def _validateInputs(self, _flnc_fa, _nfl_fa, _ccs_fofn):
        """Validate input files and return absolute expaneded paths."""
        flnc_fa, nfl_fa, ccs_fofn = _flnc_fa, _nfl_fa, _ccs_fofn
        self.add_log("Checking input files.", level=logging.INFO)
        if flnc_fa is None or nfl_fa is None:
            raise ClusterException(
                "Input non-chimeric reads files needs to be specified.")
        else:
            flnc_fa, nfl_fa = realpath(flnc_fa), realpath(nfl_fa)
            if not op.exists(flnc_fa):
                raise ClusterException("Unable to find full-length " +
                    "non-chimeric reads: {fn}".format(fn=flnc_fa))
            if not op.exists(nfl_fa):
                raise ClusterException("Unable to find non-full-length " +
                    "non-chimeric reads: {fn}".format(fn=nfl_fa))
            if ccs_fofn is not None and not op.exists(ccs_fofn):
                raise ClusterException("Unable to find FOFN file: " +
                    "{fn}".format(fn=ccs_fofn))
        return (flnc_fa, nfl_fa, ccs_fofn)

    def _validateOutputs(self, _root_dir, _out_fa):
        """Validate outputs, create root_dir if it does not exist."""
        self.add_log("Checking outputs.", level=logging.INFO)
        root_dir, out_fa = _root_dir, _out_fa
        if root_dir is None:
            self.add_log("Output directory needs to be specified.",
                         level=logging.ERROR)
        if out_fa is None:
            self.add_log("Output consensus fasta needs to be specified.",
                         level=logging.ERROR)

        root_dir = realpath(root_dir)
        out_fa = realpath(out_fa)

        if op.exists(root_dir):
            self.add_log("Output directory {d} already exists.".
                format(d=root_dir))
        else:
            self.add_log("Creating output directory {d}.".format(d=root_dir))
            os.mkdir(root_dir)
        if op.exists(out_fa):
            raise ClusterException("Consensus FASTA file {f} already exists.".
                format(f=out_fa))
        return root_dir, out_fa

    def sanity_check(self):
        """Do sanity check before stat to run."""
        errMsg = ""
        if self.ice_opts.quiver is True:
            if self.bas_fofn is None:
                errMsg = "A fofn of bas/bax.h5 files, e.g., input.fofn, " + \
                         "is required in order to polish consensus " + \
                         "isoforms using quiver."
            if self.nfl_fa is None:
                errMsg = "Non-full-length reads are required for polishing " + \
                         "consensus isoforms using quiver."
        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise ValueError(errMsg)

    @property
    def configFN(self):
        """Return configuration file of the current run."""
        return op.join(self.root_dir, "run_ice_config.txt")

    def _logConfigs(self):
        """Log configuration."""
        with open (self.configFN, 'w', 0) as f:
            f.write('pbtranscript ' + get_version() + "\n")
            f.write(str(self.ice_opts) + "\n")
            f.write(str(self.sge_opts) + "\n")

    @property
    def initPickleFN(self):
        """Return path to pickle file with initial clusters."""
        return op.join(self.root_dir, "init.uc.pickle")

    def _setProbQV(self, ccs_fofn=None, firstSplitFa=None):
        """Set probability and quality values.
        If a fofn file is specified, load QV from it. Otherwise, use
        a pre-defined probability model."""
        if ccs_fofn is None:
            self.add_log("Loading predefined probabilities model.")
            self._probqv = ProbFromModel(0.01, 0.07, 0.06)
        else:
            self.add_log("Loading probabilities and QV from {f}.".
                format(f=firstSplitFa))
            self._probqv = ProbFromQV(ccs_fofn, firstSplitFa)

    def writeSummary(self, fa, summary_fn, hq_fa=None, lq_fa=None):
        """Extract number of consensus isoforms predicted, and total
        number of bases in all consensuus isoforms from fa and write
        the two attributes to summary_fn.

        if hq_fa (polished high-quality isoforms) is not None, report
            the number of polished hq clusters
        if lq_fa (polished high-quality isoforms) is not None, report
            the number of polished hq clusters
        """
        try:
            with FastaReader(fa) as reader:
                for r in reader:
                    self.summary.numConsensusIsoforms += 1
                    self.summary.numTotalBases += len(r.sequence)
            if hq_fa is not None:
                self.summary.num_polished_hq_isoforms = 0
                with FastaReader(hq_fa) as reader:
                    for r in reader:
                        self.summary.num_polished_hq_isoforms += 1
            if lq_fa is not None:
                self.summary.num_polished_lq_isoforms = 0
                with FastaReader(lq_fa) as reader:
                    for r in reader:
                        self.summary.num_polished_lq_isoforms += 1
            self.summary.write(summary_fn)
        except ZeroDivisionError:
            errMsg = "No consensus isoforms predicted."
            self.add_log(errMsg, level=logging.ERROR)
            raise ClusterException(errMsg)

    def run(self):
        """Call ICE to cluster consensus isoforms."""
        self.add_log("Start to run cluster.", level=logging.INFO)

        # Split flnc_fa into smaller files and save files to _flnc_splitted_fas.
        self.add_log("Splitting {flnc} into ".format(flnc=self.flnc_fa) +
                     "smaller files each containing {n} reads.".format(
                     n=self.ice_opts.flnc_reads_per_split),
                     level=logging.INFO)
        self._flnc_splitted_fas = splitFasta(
            input_fasta=self.flnc_fa,
            reads_per_split=self.ice_opts.flnc_reads_per_split,
            out_dir=self.root_dir,
            out_prefix="input.split")
        self.add_log("Splitted files are: " +
                     "\n".join(self._flnc_splitted_fas),
                     level=logging.INFO)

        firstSplit = self._flnc_splitted_fas[0]
        # Set up probabbility and quality value model
        self._setProbQV(ccs_fofn=self.ccs_fofn, firstSplitFa=firstSplit)

        # Initialize cluster by clique
        # check if init.pickle already exists, if so, no need to run IceInit
        if os.path.exists(self.initPickleFN):
            self.add_log("{0} already exists. Reading to get uc.".format(self.initPickleFN), level=logging.INFO)
            with open(self.initPickleFN) as f:
                uc = cPickle.load(f)
        else:
            self.add_log("Finding maximal cliques.", level=logging.INFO)
            self.iceinit = IceInit(readsFa=firstSplit,
                      qver_get_func=self._probqv.get_smoothed,
                      ice_opts=self.ice_opts,
                      sge_opts=self.sge_opts)
            uc = self.iceinit.uc
            # Dump uc to a file
            self.add_log("Dumping initial clusters to {f}".format(
                         f=self.initPickleFN), level=logging.INFO)
            with open(self.initPickleFN, 'w') as f:
                cPickle.dump(uc, f)

        # Run IceIterative.
        self.add_log("Iteratively clustering.", level=logging.INFO)
        self.icec = IceIterative(
                fasta_filename=firstSplit,
                fasta_filenames_to_add=self._flnc_splitted_fas[1:],
                all_fasta_filename=self.flnc_fa,
                ccs_fofn=self.ccs_fofn,
                root_dir=self.root_dir,
                ice_opts=self.ice_opts,
                sge_opts=self.sge_opts,
                uc=uc,
                probQV=self._probqv)
        self.icec.run()
        clean_up_after_ICE(self.root_dir)

        # IceIterative done, write predicted (unplished) consensus isoforms
        # to an output fasta
        self.add_log("Creating a link to unpolished consensus isoforms.")
        ln(self.icec.final_consensus_fa, self.out_fa)

        # Call quiver to polish predicted consensus isoforms.
        if self.ice_opts.quiver is not True:
            self.add_log("Creating a link to cluster report.")
            ln(src=self.icec.report_fn, dst=self.report_fn)

            self.add_log("Writing a summary to {f}".format(f=self.summary_fn),
                         level=logging.INFO)
            self.writeSummary(fa=self.out_fa, summary_fn=self.summary_fn)
        else:  # self.ice_opts.quiver is True
            #TODO review code
            self.pol = Polish(root_dir=self.root_dir,
                         nfl_fa=self.nfl_fa,
                         bas_fofn=self.bas_fofn,
                         ccs_fofn=self.ccs_fofn,
                         hq_isoforms_fa=self.hq_isoforms_fa,
                         hq_isoforms_fq=self.hq_isoforms_fq,
                         lq_isoforms_fa=self.lq_isoforms_fa,
                         lq_isoforms_fq=self.lq_isoforms_fq,
                         ice_opts=self.ice_opts,
                         sge_opts=self.sge_opts)
            self.pol.run()

            # cluster report
            self.add_log("Creating a link to cluster report.")
            ln(src=self.pol.iceq.report_fn, dst=self.report_fn)

            # Write a summary.
            self.add_log("Writing a summary to {f}".format(f=self.summary_fn),
                         level=logging.INFO)
            self.writeSummary(fa=self.out_fa, summary_fn=self.summary_fn,
                              hq_fa=self.pol.icepq.quivered_good_fa,
                              lq_fa=self.pol.icepq.quivered_bad_fa)

        # Create log file.
        self.close_log()
        return 0


