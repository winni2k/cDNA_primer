#!/usr/bin/env python
"""
For each input fasta file, map its reads to consensus isoforms
in ref_fasta and create a pickle file (e.g. *.partial_uc.pickle).
When all pickle files are created, merge them and dump to a big pickle.
"""

from collections import defaultdict
from cPickle import dump, load
import os.path as op
import os
import logging
import sys
import shutil
import time
from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.PBTranscriptOptions import add_sge_options, \
        add_fofn_arguments
from pbtools.pbtranscript.Utils import realpath, mkdir
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ClusterOptions import SgeOptions
from pbcore.util.Process import backticks

#log = logging.getLogger(__name__)

class IceAllPartials(IceFiles):
    """Include non-full-length sequences and build partial uc.
    """
    def __init__(self, root_dir, fasta_filenames, ref_fasta,
            out_pickle, sge_opts, sa_file=None, ccs_fofn=None):
        """
        fasta_filenames --- a list of splitted nfl fasta files.

        ref_fasta --- (unpolished) consensus isoforms

        out_pickle --- a pickle file with all nfl fasta reads

        ccs_fofn --- should be reads_of_insert.fofn or None

        root_dir --- ICE root output directory

        sge_opts --- params for SGE environment, including
            use_sge    : use SGE or not
            max_sge_jobs: maximum number of gcon jobs submitted
            unique_id  : unique qsub job id, important that this
                        DOES NOT CONFLICT!
            blasr_nproc: blasr -nproc param, number of threads per cpu.
            gcon_nproc : number of gcon that can run at the same time
        """
        self.prog_name = "IceAllPartials"
        IceFiles.__init__(self, prog_name=self.prog_name, root_dir=root_dir)

        self.fasta_filenames, self.ref_fasta, self.ccs_fofn, self.sa_file = \
            self._validateInputs(fasta_filenames=fasta_filenames,
                                 ref_fasta=ref_fasta,
                                 ccs_fofn=ccs_fofn,
                                 sa_file=sa_file)

        self.out_pickle = out_pickle

        self.sge_opts = sge_opts

        self.add_log("Making dir for mapping noFL reads: " + self.nfl_dir)
        mkdir(self.nfl_dir)

        self.add_log("input fasta files are: " +
                    ", ".join(self.fasta_filenames))
        self.add_log("temp pickle files are: " +
                    ", ".join(self.pickle_filenames))
        self.add_log("out pickle file is: " + self.out_pickle)

    def _validateInputs(self, fasta_filenames, ref_fasta, ccs_fofn, sa_file):
        """Validate input files."""
        for f in fasta_filenames:
            if not op.exists(f):
                raise IOError("Input fasta {f} does not exist.".format(f=f))
        if ref_fasta is None or not op.exists(ref_fasta):
            raise IOError("Reference {r} does not exist.".format(r=ref_fasta))
        if ccs_fofn is not None and not op.exists(ccs_fofn):
            raise IOError("ccs_fofn file {ccs_fofn} does not exist.".format(
                ccs_fofn=ccs_fofn))
        if sa_file is not None and not op.exists(sa_file):
            raise IOError("suffix array {s} does not exist.".format(s=sa_file))
        return ([realpath(f) for f in fasta_filenames],
                realpath(ref_fasta),
                realpath(ccs_fofn),
                realpath(sa_file) if sa_file is not None else None)

    @property
    def pickle_filenames(self):
        """pickle files for each fasta file."""
        return [op.join(self.nfl_dir, op.basename(f) + ".partial_uc.pickle")
                for f in self.fasta_filenames]

    @property
    def done_filenames(self):
        """done files to indicate that pickles are done."""
        return [op.join(self.nfl_dir, op.basename(f) + ".DONE")
                for f in self.pickle_filenames]

    @property
    def script_filenames(self):
        """scripts to generate pickles from fasta files."""
        return [op.join(self.script_dir, op.basename(f) + ".partial_uc.sh")
                for f in self.fasta_filenames]

    def createPickles(self):
        """For each file in fasta_filenames, call ice_partial.py to build
        clusters and to save results to a pickle file. When all pickles
        are done, union all pickles.
        """
        self.add_log("Mapping non-full-length reads to consensus isoforms.")
        self.add_log("Creating pickles...", level=logging.INFO)

        for idx, fa in enumerate(self.fasta_filenames):
            # for each splitted non-full-length reads fasta file, build #
            # partial_uc.pickle
            self.add_log("Creating a pickle for {f}".format(f=fa))
            cmd = "ice_partial.py {i} ".format(i=fa) + \
                  "{r} ".format(r=self.ref_fasta) + \
                  "{o} ".format(o=self.pickle_filenames[idx]) + \
                  "--blasr_nproc={n} ".format(n=self.sge_opts.blasr_nproc) + \
                  "--done={d} ".format(d=self.done_filenames[idx])
            if self.ccs_fofn is not None:
                cmd += "--ccs_fofn={f} ".format(f=self.ccs_fofn)
            if self.sa_file is not None:
                cmd += "--sa={sa} ".format(sa=self.sa_file)

            self.add_log("Writing command to script {fsh}".
                         format(fsh=self.script_filenames[idx]))
            with open(self.script_filenames[idx], 'w') as fsh:
                fsh.write(cmd + "\n")

            # determine elog & olog
            partial_log_fn = op.join(self.log_dir,
                                     'IcePartial.{idx}'.format(idx=idx))
            elog = partial_log_fn + ".elog"
            olog = partial_log_fn + ".olog"
            jid = "ice_partial_" + op.basename(fa)

            qsub_cmd = "qsub " + \
                       "-pe smp {n} ".format(n=self.sge_opts.blasr_nproc) + \
                       "-cwd -S /bin/bash -V " + \
                       "-e {elog} ".format(elog=elog) + \
                       "-o {olog} ".format(olog=olog) + \
                       "-N {jid} ".format(jid=jid) + \
                       "{sh}".format(sh=self.script_filenames[idx])

            if self.sge_opts.use_sge is True:
                self.add_log("Submitting CMD: {qcmd}".format(qcmd=qsub_cmd))
                _out, _code, _msg = backticks(qsub_cmd)
            #elif self.sge_opts.useSMRTPortal is True:
            #    pass
            else:
                cmd += " 1>{olog} 2>{elog}".format(olog=olog, elog=elog)
                self.add_log("Submitting CMD: {cmd}".format(cmd=cmd))
                _out, _code, _msg = backticks(cmd)
                if _code != 0:
                    raise RuntimeError("CMD failed: {cmd}\n{msg}\n".format(
                        cmd=cmd, msg=str(_msg)))

    def waitForPickles(self, pickle_filenames, done_filenames):
        """Wait for *.pickle and *.pickle.DONE to be created."""
        self.add_log("Waiting for pickles {ps} to be created.".
                     format(ps=", ".join(pickle_filenames)),
                     level=logging.INFO)

        stop = False
        sleep_time = 10
        while stop is not True:
            stop = all(op.exists(p) for p in pickle_filenames) and \
                   all(op.exists(d) for d in done_filenames)
            sleep_time = min(600, sleep_time * 2)
            time.sleep(sleep_time)
            self.add_log("Waiting for pickles to be created: {ps}".format(ps=
                        ", ".join([p for p in pickle_filenames
                                   if op.exists(p)])))

    def combinePickles(self, pickle_filenames, out_pickle):
        """Combine all *.pickle files to one and dump to self.out_pickle."""
        self.add_log("Combining pickles: {ps} to a big pickle {p}".format(
                     ps=", ".join(pickle_filenames), p=out_pickle),
                     level=logging.INFO)
        if len(pickle_filenames) == 1:
            src = pickle_filenames[0]
            dst = out_pickle

            if (realpath(src) != realpath(dst)):
                self.add_log("Copying {src} to {dst}.".format(src=src, dst=dst))
                shutil.copyfile(src, dst)
            else:
                self.add_log("{dst} has been created, no need to merge.".
                    format(dst=out_pickle))
        else:
            # Combine all partial outputs
            self.add_log("Merging all pickles.")
            partial_uc = defaultdict(lambda: [])
            nohit = set()
            for pf in pickle_filenames:
                self.add_log("Merging {pf}.".format(pf=pf))
                a = load(open(pf))
                nohit.update(a['nohit'])
                for k, v in a['partial_uc'].iteritems():
                    partial_uc[k] += v

            self.add_log("Dumping all to {f}".format(f=out_pickle))
            # Dump to one file
            partial_uc = dict(partial_uc)
            with open(out_pickle, 'w') as f:
                dump({'nohit':nohit, 'partial_uc':partial_uc}, f)

    def run(self):
        """Assigning nfl reads to consensus isoforms and merge."""
        # Call ice_partial.py to create a pickle for each splitted nfl fasta
        self.createPickles()
        # Wait for pickles to be created, if SGE is used.
        self.waitForPickles(pickle_filenames=self.pickle_filenames,
                            done_filenames=self.done_filenames)
        # Combine all pickles to a big pickle file: nfl_all_pickle_fn.
        self.combinePickles(pickle_filenames=self.pickle_filenames,
                            out_pickle=self.nfl_all_pickle_fn)
        # Create symbolic link if necessary
        if realpath(self.nfl_all_pickle_fn) != realpath(self.out_pickle):
            self.add_log("Creating a symbolic link for {f}".format(
                f=self.out_pickle), level=logging.INFO)
            if op.exists(self.out_pickle):
                os.remove(self.out_pickle)
            os.symlink(self.nfl_all_pickle_fn, self.out_pickle)

        # Close log
        self.close_log()


def set_parser(parser):
    """Get argument parser."""
    parser.add_argument("fasta_filenames",
                        type=str,
                        help="comma delimited fasta files of " +
                             "splitted non-full-length reads")
    parser.add_argument("ref_fasta",
                        type=str,
                        help="Reference fasta file, most likely " +
                             "ref_consensus.fa from ICE output")
    parser.add_argument("out_pickle",
                        type=str,
                        help="Output pickle file.")
    parser.add_argument("--root_dir",
                        dest="root_dir",
                        default="",
                        help="A directory for saving intermediate files.")

    parser = add_sge_options(parser, blasr_nproc=True)
    parser = add_fofn_arguments(parser, ccs_fofn=True)

    parser.add_argument("--sa",
                        dest="sa_file",
                        default=None,
                        help="saffix array of ref_fasta")
    return parser


from pbcore.util.ToolRunner import PBToolRunner

class IceAllPartialsRunner(PBToolRunner):
    """IceAllPartials Runner"""
    def __init__(self):
        desc = "Assign all non-full-length reads to unpolished isoform clusters"
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
        sge_opts = SgeOptions(os.getpid(),
                blasr_nproc=args.blasr_nproc,
                use_sge=args.use_sge,
                max_sge_jobs=args.max_sge_jobs)

        icep = IceAllPartials(
                root_dir=args.root_dir,
                fasta_filenames=args.fasta_filenames.split(','),
                ref_fasta=args.ref_fasta,
                out_pickle=args.out_pickle,
                sge_opts=sge_opts,
                sa_file=args.sa_file,
                ccs_fofn=args.ccs_fofn)

        icep.run()
        return 0


def main():
    """Main function"""
    runner = IceAllPartialsRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
