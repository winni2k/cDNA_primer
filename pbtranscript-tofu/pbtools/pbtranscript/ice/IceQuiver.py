#!/usr/bin/env python

"""Call quiver to polish consensus isoforms created by ICE."""

import os, sys
import os.path as op
import logging
from cPickle import load
from collections import defaultdict
from pbcore.util.Process import backticks
from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.PBTranscriptOptions import  add_fofn_arguments, \
        add_sge_options
from pbtools.pbtranscript.Utils import mkdir
from pbtools.pbtranscript.ice.IceUtils import get_the_only_fasta_record, \
        get_files_from_fofn, is_blank_sam, concat_sam, \
        blasr_sam_for_quiver, write_in_raw_fasta
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.io.FastaRandomReader import MetaSubreadFastaReader


class IceQuiver(IceFiles):
    """Ice Quiver."""
    def __init__(self, root_dir, bas_fofn, fasta_fofn, sge_opts):
        #Initialize super class IceFiles.
        super(IceQuiver, self).__init__(prog_name="IceQuiver",
            root_dir=root_dir, bas_fofn=bas_fofn, fasta_fofn=fasta_fofn)

        self.sge_opts = sge_opts
        self._validateInputs()

    def _validateInputs(self):
        """Validate input fofns, and root_dir, log_dir, tmp_dir"""
        self.add_log("Validating inputs.")
        errMsg = ""

        if not op.exists(self.log_dir) or not op.isdir(self.log_dir):
            errMsg = "Log dir {l} is not an existing directory !".\
                format(l=self.log_dir)
        if not op.exists(self.bas_fofn):
            errMsg = "Failed to find fofn {f} for bas/bax.h5 files.".\
                format(f=self.bas_fofn)
        if not op.exists(self.fasta_fofn):
            errMsg = "Failed to find fasta fofn {f}!".format(f=self.fasta_fofn)
        if not op.exists(self.nfl_all_pickle_fn):
            #"output/map_noFL/noFL.ALL.partial_uc.pickle"):
            errMsg = "Failed to find {f}!".format(f=self.nfl_all_pickle_fn)

        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise ValueError(errMsg)

    def _quivered_bin_prefix(self, first, last):
        """Return $quivered_dir/c{first}to{last}"""
        return self.quivered_dir + "/c{first}to{last}".format(
                first=first, last=last)

    def sam_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.sam"""
        return self._quivered_bin_prefix(first, last) + ".sam"

    def ref_fa_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.ref.fa
        this is reference fasta for quiver to use as input.
        """
        #
        return self._quivered_bin_prefix(first, last) + ".ref.fa"

    def cmph5_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.cmp.h5"""
        return self._quivered_bin_prefix(first, last) + ".cmp.h5"

    def fq_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.quivered.fq
        this is quivered fq output. Whenever this is changed, change
        IcePostQuiver accordingly.
        """
        return self._quivered_bin_prefix(first, last) + ".quivered.fq"

    def script_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.sh"""
        return self._quivered_bin_prefix(first, last) + ".sh"

    def setup_quiver_for_batch(self, cids, refs, quiver_nproc=2,
            return_script=True):
        """
        NOTE: (1) skip clusters if identical sequences already exists
                  in another cluster (rare, but happens)
              (2) skip clusters if the alignment is empty (also rare,
                  but happens)

        if return_script is True, return a job script e.g., quivered/c{}to{}.sh,
        otherwise, return a list of cmds.
        """
        # concat the sam files
        first, last = cids[0], cids[-1]
        #prefix = self._quivered_bin_prefix(first=cids[0], last=cids[-1])
        bin_sam_file = self.sam_of_quivered_bin(first, last)
        bin_ref_fa = self.ref_fa_of_quivered_bin(first, last)
        bin_cmph5 = self.cmph5_of_quivered_bin(first, last)
        bin_fq = self.fq_of_quivered_bin(first, last)

        valid_sam_files = []
        valid_cids = []
        seqs_seen = {}
        for cid in cids:
            fname = self.sam_of_cluster(cid)
            if not is_blank_sam(fname):
                seq = get_the_only_fasta_record(refs[cid]).sequence
                if seq not in seqs_seen:
                    valid_sam_files.append(fname)
                    valid_cids.append(cid)
                    seqs_seen[seq] = cid
                else:
                    self.add_log(
                        "ignoring {0} because identical sequence!".format(cid))
            else:
                self.add_log(
                    "ignoring {0} because no alignments!".format(cid))
        concat_sam(valid_sam_files, bin_sam_file)

        # concat the reference file
        cmd = "cat " + " ".join(refs[cid] for cid in valid_cids) + \
              " > {ref}".format(ref=bin_ref_fa)
        _out, _code, _msg = backticks(cmd)
        if _code != 0:
            errMsg = "Unable to concatenate reference files between " + \
                "{first} and {last}.\n".format(first=first, last=last) + _msg
            self.add_log(errMsg, level=logging.ERROR)
            raise RuntimeError(errMsg)

        # write the sh script for the conversion, loadPulses, and quiver
        cmds = []
        cmds.append("samtoh5 {sam} {ref} {cmph5} -smrtTitle".format(
            sam=bin_sam_file, ref=bin_ref_fa, cmph5=bin_cmph5))
        cmds.append("gzip {sam}".format(sam=bin_sam_file))
        metrics = ["QualityValue", "InsertionQV", "MergeQV", "DeletionQV",
                   "DeletionTag", "SubstitutionTag", "SubstitutionQV"]
        cmds.append("loadPulses {bas_fofn} ".format(bas_fofn=self.bas_fofn) +
                    "{cmph5} ".format(cmph5=bin_cmph5) +
                    "-byread -metrics " + ",".join(metrics))
        cmds.append("cmph5tools.py sort {cmph5}".format(cmph5=bin_cmph5))
        cmds.append("samtools faidx {ref}".format(ref=bin_ref_fa))
        cmds.append("quiver {cmph5} ".format(cmph5=bin_cmph5) +
                    "-v -j{n} ".format(n=quiver_nproc) +
                    "-r {ref} ".format(ref=bin_ref_fa) +
                    "-o {fq}".format(fq=bin_fq))

        bin_sh = self.script_of_quivered_bin(first, last)

        with open(bin_sh, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("\n".join(cmds))

        if return_script is True:
            return f.name
        else:
            return cmds

    def submit_quiver_jobs(self, d, uc, partial_uc, refs, keys, start, end,
                           submitted, todo,
                           use_sge, max_sge_jobs, quiver_nproc):
        """Call quiver to polish consensus.
        (1) for each cluster k, obtain unrolled sequences of all reads (zmws)
            belonging to this cluster, and save in raw_fa_of_cluster(k)
        (2) for each cluster k, call blasr to align raw_f_of_cluster to
            consensus sequence of the cluster and create sam_of_cluster.

        (3) Put every 100 clusters into one big bin, and then
            merge all sam_of_cluster files to sam_of_quivered_bin
        (4) Prepare commands including
                samtoh5, loadPulses, cmph5tools.py ...
            in order to convert sam_of_quivered_bin to cmph5_of_quivered_bin.
                * Either write these command to script_of_quivered_bin and qsub
                  all jobs later when scripts of all quivered bins are done,
                * Or execute the commands immediately.
        """
        for i in xrange(start, end, 100):
            for k in keys[i: min(end, i+100)]:
                #os.chdir(op.join('./tmp', str(k/10000), 'c'+str(k)))
                raw_fa = self.raw_fa_of_cluster(k)

                # write_in_raw_fa return movies of reads in partial_uc
                # logging.debug("uc[k]={0}".format(uc[k]))
                # logging.debug("partial_uc[k]={0}".format(uc[k]))
                write_in_raw_fasta(input_fasta_d=d,
                    in_seqids=uc[k] + partial_uc[k],
                    out_fa=raw_fa,
                    ignore_keyerror=True)

                #TODO: use multi-processing pool, reduce nproc
                blasr_sam_for_quiver(
                    input_fasta=raw_fa,
                    ref_fasta=refs[k],
                    out_sam_filename=self.sam_of_cluster(k),
                    run_cmd=True)

            fname = self.setup_quiver_for_batch(cids=keys[i: min(end, i+100)],
                        refs=refs, quiver_nproc=quiver_nproc,
                        return_script=True)
            todo.append(fname)

            if use_sge is not True or \
               max_sge_jobs == 0: # don't use SGE
                for job in todo:
                    elog = op.join(self.quivered_log_dir,
                                   op.basename(job) + ".elog")
                    olog = op.join(self.quivered_log_dir,
                                   op.basename(job) + ".olog")
                    msg = "Running quiver job locally: {j} ".format(j=job) + \
                          "1>{olog} 2>{elog}".format(olog=olog, elog=elog)
                    self.add_log(msg)
                    cmd = "bash " + job + " 1>{olog} 2>{elog}".\
                          format(olog=olog, elog=elog)
                    _out, _code, _msg = backticks(cmd)
                    if _code != 0:
                        errMsg = "Failed to run quiver {j}".format(j=job) + _msg
                        self.add_log(errMsg, level=logging.ERROR)
                        raise RuntimeError(errMsg)
                    submitted.append(("local", job))
                todo = []
            else:
                while len(todo) > 0:
                    n = min(max_sge_jobs, len(todo))
                    for job in todo[:n]:
                        # ex: Your job 8613116 ("c20to70.sh") has been submitted
                        elog = op.join(self.quivered_log_dir,
                                       op.basename(job) + ".elog")
                        olog = op.join(self.quivered_log_dir,
                                       op.basename(job) + ".olog")
                        qsub_cmd = "qsub " + \
                                   "-pe smp {n} ".format(n=quiver_nproc) + \
                                   "-cwd -S /bin/bash -V " + \
                                   "-e {elog} ".format(elog=elog) + \
                                   "-o {olog} ".format(olog=olog) + \
                                   "{job}".format(job=job)
                        msg = "Submitting CMD: {cmd}.\n".format(cmd=qsub_cmd)
                        self.add_log(msg)
                        _out, _code, _msg = backticks(qsub_cmd)
                        if _code != 0:
                            errMsg = "Failed to submit CMD {cmd}.".format(
                                    cmd=qsub_cmd)
                            self.add_log(errMsg, level=logging.ERROR)
                            raise RuntimeError(errMsg)

                        job_id = str(_out).split()[2]
                        submitted.append((job_id, job))
                        todo.remove(job)
                    # end of for job in todo[:n]
                # end of while len(todo) > 0
            # end of else (use sge)
        # end of for i in xrange(start, end, 100):

    @property
    def report_fn(self):
        """Return a csv report with cluster_id, read_id, read_type."""
        return op.join(self.out_dir, "cluster_report.FL_nonFL.csv")

    def write_report(self, uc, partial_uc, report_fn):
        """
        Write a CSV report to report_fn, each line contains three columns:
            cluster_id, read_id and read_type
        """
        with open(report_fn, 'w') as f:
            f.write("cluster_id\tread_id\tread_type\n")
            for c in uc.keys():
                for r in uc[c]:
                    f.write("c{c}\t{r}\tFL\n".format(r=r, c=c))
                for r in partial_uc[c]:
                    f.write("c{c}\t{r}\tNonFL\n".format(r=r, c=c))

    def run(self):
        """Run quiver for ICE."""
        # Create directories: root_dir/quivered and root_dir/log_dir/quivered
        mkdir(self.quivered_dir)
        mkdir(self.quivered_log_dir)

        files = get_files_from_fofn(self.fasta_fofn)
        msg = "Indexing {0} fasta files, please wait.".format(len(files))
        self.add_log(msg)

        d = MetaSubreadFastaReader(files)
        self.add_log("Fasta files indexing done.")

        self.add_log("Loading uc from {f}.".format(f=self.final_pickle_fn))
        a = load(open(self.final_pickle_fn))
        uc = a['uc']
        refs = a['refs']

        self.add_log("Loading partial uc from {f}.".
                     format(f=self.nfl_all_pickle_fn))
        partial_uc = load(open(self.nfl_all_pickle_fn))['partial_uc']
        partial_uc2 = defaultdict(lambda: [])
        partial_uc2.update(partial_uc)

        # Write report to quivered/cluster_report.FL_nonFL.csv
        self.add_log("Writing a csv report of cluster -> FL/NonFL reads to {f}".
                     format(f=self.report_fn), level=logging.INFO)
        self.write_report(uc=uc, partial_uc=partial_uc2,
                          report_fn=self.report_fn)

        good = [x for x in uc] #[x for x in uc if len(uc[x]) > 1 or len(partial_uc2[x]) >= 10]
        keys = sorted(good)  # sort good keys (cluster ids)

        start = 0
        end = len(keys)

        submitted = []  # submitted jobs
        todo = []       # to-do jobs

        self.submit_quiver_jobs(d=d, uc=uc, partial_uc=partial_uc2,
            refs=refs, keys=keys, start=start, end=end,
            submitted=submitted, todo=todo,
            use_sge=self.sge_opts.use_sge,
            max_sge_jobs=self.sge_opts.max_sge_jobs,
            quiver_nproc=self.sge_opts.quiver_nproc)

        with open(self.submitted_quiver_jobs_log, 'w') as f:
            f.write("\n".join(str(x[0]) + '\t' + str(x[1]) for x in submitted))

        self.close_log()
        return 0


def set_parser(parser):
    """Set argument parser."""
    parser.add_argument("root_dir",
                        type=str,
                        default=None,
                        help="Root directory for ICE, " + \
                             "e.g., SMRTPortal data/clusterOutDir")
    parser = add_fofn_arguments(parser, bas_fofn=True, fasta_fofn=True)
    parser = add_sge_options(parser, quiver_nproc=True)


from pbcore.util.ToolRunner import PBToolRunner
class IceQuiverRunner(PBToolRunner):
    """IceQuiver runner."""
    def __init__(self):
        desc = "Call quiver to polish consensus isoforms."
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
        from pbtools.pbtranscript.ClusterOptions import SgeOptions
        sge_opts = SgeOptions(os.getpid(),
                              use_sge=args.use_sge,
                              max_sge_jobs=args.max_sge_jobs,
                              quiver_nproc=args.quiver_nproc)
        try:
            obj = IceQuiver(root_dir=args.root_dir,
                            bas_fofn=args.bas_fofn,
                            fasta_fofn=args.fasta_fofn,
                            sge_opts=sge_opts)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            import traceback
            traceback.print_exc()
            return 1
        return 0

def main():
    """Main function"""
    runner = IceQuiverRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
