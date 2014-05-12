#!/usr/bin/env python

"""Polish clustered isoforms using quiver."""

import sys
import os.path as op
import logging
from pbtools.pbtranscript.io.FastaSplitter import splitFasta
from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.Utils import realpath
from pbtools.pbtranscript.ice.IceUtils import build_sa, \
        convert_fofn_to_fasta
from pbtools.pbtranscript.PBTranscriptOptions import add_fofn_arguments, \
        add_polished_isoforms_arguments, add_sge_options
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ice.IceAllPartials import IceAllPartials
from pbtools.pbtranscript.ice.IceQuiver import IceQuiver
from pbtools.pbtranscript.ice.IcePostQuiver import IcePostQuiver


class Polish(IceFiles):
    """Polish isoforms clusters using Quiver."""
    def __init__(self, root_dir, nfl_fa, bas_fofn, ccs_fofn,
                 ice_opts, sge_opts, hq_isoforms_fa=None, hq_isoforms_fq=None,
                 lq_isoforms_fa=None, lq_isoforms_fq=None, fasta_fofn=None):
        """
        root_dir --- IceFiles.root_dir, usually data/clusterOutDir
        nfl_fa    --- non-full-length reads in fasta, e.g., isoseq_nfl.fa
        bas_fofn --- e.g. input.fofn of bas|bax.h5 files
        ccs_fofn --- e.g. reads_of_insert.fofn of ccs files.
        hq_isoforms_fa|fq  --- polished, hiqh quality consensus isoforms in fasta|q
        lq_isoforms_fa|fq  --- polished, low quality consensus isoforms in fasta|q
        """
        IceFiles.__init__(self, prog_name="IcePolish", root_dir=root_dir,
                          bas_fofn=bas_fofn, ccs_fofn=ccs_fofn,
                          fasta_fofn=fasta_fofn)
        self.nfl_fa = realpath(nfl_fa)
        self.hq_isoforms_fa = hq_isoforms_fa
        self.hq_isoforms_fq = hq_isoforms_fq
        self.lq_isoforms_fa = lq_isoforms_fa
        self.lq_isoforms_fq = lq_isoforms_fq
        self.ice_opts = ice_opts
        self.sge_opts = sge_opts

        self.icep = None   # IceAllPartials.
        self.iceq = None   # IceQuiver
        self.icepq = None  # IcePostQuiver
        self._nfl_splitted_fas = None

        self._validate_inputs()

    def _validate_inputs(self):
        """
        Validate input directories: root_dir, and
        files: nfl_fa, bas_fofn, ccs_fofn.
        """
        self.add_log("Validating inputs.")
        errMsg = ""
        if not op.exists(self.root_dir):
            errMsg = "Root dir {d} is not an existing directory!".\
                format(d=self.root_dir)
        if not op.exists(self.nfl_fa):
            errMsg = "Failed to find non-full-length reads {f}!".\
                format(f=self.nfl_fa)
        if not op.exists(self.bas_fofn):
            errMsg = "Failed to find bas fofn {f}!".format(f=self.bas_fofn)
        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise ValueError(errMsg)

    def get_sa_file(self):
        """Return saffix array of final_consensus_fa."""
        ret = True
        if not op.exists(self.final_consensus_sa):
            ret = build_sa(self.final_consensus_fa, self.final_consensus_sa)

        if ret is True:
            self.add_log("Successfully created a suffix array for {f} at {sa}.".
                         format(f=self.final_consensus_fa,
                                sa=self.final_consensus_sa))
        else:
            self.add_log("Failed to generate suffix array for {f}".format(
                f=self.final_consensus_fa), level=logging.WARNING)
        if op.exists(self.final_consensus_sa):
            return self.final_consensus_sa
        else:
            return None

    def run(self):
        """
        First, split non-full-length (nfl) fasta files into smaller
        chunks, assign nfl reads in each splitted fasta file
        into unpolished isoform clusters and then merge all pickles
        into self.nfl_all_pickle_fn.
        Second, bin every 100 clusters, for each bin, call blasr,
        samto5h, loadPulses, cmph5tools to create cmp.h5 files and
        call quiver to polish each isoforms within each bin.
        Finally, pick up good isoform clusters whose QV errors is less
        than a threshold.
        Save all high quality isoforms to hq_isoforms_fa|fq if they are not None
        Save all low quality isoforms to lq_isoforms_fa|fq if they are not None
        """
        # Create final.consensus.fa.sa
        self.add_log("Generating suffix array for {f}".format(
                     f=self.final_consensus_sa), level=logging.INFO)
        sa_file = self.get_sa_file()

        # Create input.fasta.fofn from bas_fofn
        self.add_log("Creating fasta fofn from bas/bax.h5 fofn",
                     level=logging.INFO)
        if self.fasta_fofn is None:
            self.fasta_fofn = op.join(self.nfl_dir, "input.fasta.fofn")
        self.add_log("bas fofn={f}".format(f=self.bas_fofn))
        self.add_log("fasta fofn={f}".format(f=self.fasta_fofn))
        convert_fofn_to_fasta(fofn_filename=self.bas_fofn,
                              out_filename=self.fasta_fofn,
                              fasta_out_dir=self.nfl_dir)

        # Split non-full-length reads into smaller fasta files
        # and save files to root_dir/nfl_00.fa, ..., .
        self.add_log("Splitting {nfl} into ".format(nfl=self.nfl_fa) +
                     "smaller files each containing {n} reads.".format(
                     n=self.ice_opts.nfl_reads_per_split),
                     level=logging.INFO)
        self._nfl_splitted_fas = splitFasta(input_fasta=self.nfl_fa,
            reads_per_split=self.ice_opts.nfl_reads_per_split,
            out_dir=self.nfl_dir,
            out_prefix="input.split")
        msg = "Splitted files are: " + "\n".join(self._nfl_splitted_fas)
        self.add_log(msg, level=logging.INFO)

        # Process nfl reads in each splitted fasta.
        self.add_log("IceAllPartials initiated.", level=logging.INFO)
        sa_file = self.final_consensus_sa \
                  if op.exists(self.final_consensus_fa) else None
        self.icep = IceAllPartials(
                root_dir=self.root_dir,
                fasta_filenames=self._nfl_splitted_fas,
                ref_fasta=self.final_consensus_fa,
                out_pickle=self.nfl_all_pickle_fn,
                sge_opts=self.sge_opts,
                sa_file=sa_file,
                ccs_fofn=self.ccs_fofn)
        self.icep.run()
        self.add_log("IceAllPartials completed.", level=logging.INFO)

        self.add_log("IceQuiver initiated.", level=logging.INFO)
        self.iceq = IceQuiver(root_dir=self.root_dir,
                              bas_fofn=self.bas_fofn,
                              fasta_fofn=self.fasta_fofn,
                              sge_opts=self.sge_opts)
        self.iceq.run()
        self.add_log("IceQuiver finished.", level=logging.INFO)

        self.add_log("IcePostQuiver initiated.", level=logging.INFO)
        self.icepq = IcePostQuiver(root_dir=self.root_dir,
                                   hq_isoforms_fa=self.hq_isoforms_fa,
                                   hq_isoforms_fq=self.hq_isoforms_fq,
                                   lq_isoforms_fa=self.lq_isoforms_fa,
                                   lq_isoforms_fq=self.lq_isoforms_fq,
                                   use_sge=self.sge_opts.use_sge,
                                   quit_if_not_done=False)
        self.icepq.run()
        self.add_log("IcePostQuiver finished.", level=logging.INFO)


def set_parser(parser):
    """Set up argument parser."""
    parser.add_argument("root_dir",
                        type=str,
                        default=None,
                        help="Root directory for ICE, " + \
                             "e.g., SMRTPortal data/clusterOutDir")
    parser = add_nfl_fa_arguments(parser, positional=True)

    parser = add_fofn_arguments(parser, ccs_fofn=True, bas_fofn=True)

    parser = add_polished_isoforms_arguments(parser)

    parser = add_sge_options(parser, quiver_nproc=True)


from pbcore.util.ToolRunner import PBToolRunner


class PolishRunner(PBToolRunner):
    """Polish Runner"""
    def __init__(self):
        desc = "Call quiver to polish consensus isoforms with " + \
               "non-full-length reads."
        PBToolRunner.__init__(self, desc)
        set_parser(self.parser)

    def getVersion(self):
        """Get version string."""
        return get_version()

    def run(self):
        """Run"""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))
        args = self.args
        from pbtools.pbtranscript.ClusterOptions import IceOptions, SgeOptions
        sge_opts = SgeOptions(unique_id=args.unique_id,
                              use_sge=args.use_sge,
                              max_sge_jobs=args.max_sge_jobs,
                              quiver_nproc=args.quiver_nproc)

        try:
            obj = Polish(root_dir=args.root_dir,
                         nfl_fa=args.nfl_fa,
                         bas_fofn=args.bas_fofn,
                         ccs_fofn=args.ccs_fofn,
                         hq_isoforms_fa=args.hq_isoforms_fa,
                         hq_isoforms_fq=args.hq_isoforms_fq,
                         lq_isoforms_fa=args.lq_isoforms_fa,
                         lq_isoforms_fq=args.lq_isoforms_fq,
                         sge_opts=sge_opts,
                         ice_opts=IceOptions())
            obj.run()
        except Exception as e:
            logging.error(str(e))
            import traceback
            traceback.print_exc()
            return 1
        return 0


def main():
    """Main function."""
    runner = PolishRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
