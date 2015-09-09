#!/usr/bin/env python
"""
Define files which are used or created in the ICE algorithm,
including temporary results, output files, scripts and logs.
"""
import os
import os.path as op
import time
import logging
from multiprocessing import Process
from pbcore.util.Process import backticks
from pbcore.io import FastaReader
from pbtools.pbtranscript.Utils import real_ppath, now_str, mkdir
from pbtools.pbtranscript.io.Summary import ClusterSummary


class IceFiles(object):

    """Define directories and files used by the ICE algorithm."""

    def __init__(self, prog_name, root_dir,
                 bas_fofn=None, ccs_fofn=None, fasta_fofn=None,
                 no_log_f=False):
        """
        prog_name --- name of a sub-class
        root_dir --- root directory of the whole project. There will be
                     sub-directories under it, including:
                     tmp/ --- 0/  c0, c1, ..., c9999
                          --- 1/  c10000, c10001, ..., c19999
                          ...
                          each c? folder contains data for a cluster id=c?
                     script/
                          --- 0/  gcon_job_?.sh, gcon jobs in the first iteration
                          --- 1/  gcon_job_?.sh, gcon jobs in the second iteration
                          ...
                     log/
                          --- ICE.log   Log of the ICE algorithm
                          --- 0/  log for jobs in the first iteration
                          ...
                     output/   output files go here.
        bas_fofn --- input.fofn which contains movie.bas|bax.h5 files.
        ccs_fofn --- a fofn contains movie.ccs.h5 files.
        fasta_fofn --- a fofn contains movie.bax.h5.fasta files.
                     script/
        no_log_f --- DON'T write log to a log file.
        """
        self.prog_name = str(prog_name)
        self.root_dir = real_ppath(root_dir)

        self.bas_fofn = real_ppath(bas_fofn)
        self.ccs_fofn = real_ppath(ccs_fofn)
        self.fasta_fofn = real_ppath(fasta_fofn)

        mkdir(self.root_dir)
        mkdir(self.tmp_dir)
        mkdir(self.log_dir)
        mkdir(self.script_dir)
        mkdir(self.out_dir)

        self.no_log_f = no_log_f
        if not no_log_f:
            self.log_f = open(self.log_fn, 'w', 0)
            self.add_log(msg="Initializing {p}.".format(p=self.prog_name))

    @property
    def tmp_dir(self):
        """Return $root_dir/tmp."""
        return op.join(self.root_dir, "tmp")

    @property
    def log_dir(self):
        """Return $root_dir/log."""
        return op.join(self.root_dir, "log")

    @property
    def log_fn(self):
        """Return $log_dir/$prog_name.log"""
        return op.join(self.log_dir, self.prog_name + ".log")

    @property
    def out_dir(self):
        """Return $root_dir/output."""
        return op.join(self.root_dir, "output")

    @property
    def script_dir(self):
        """Return $root_dir/scripts."""
        return op.join(self.root_dir, "scripts")

    @property
    def nfl_dir(self):
        """Return $root_dir/output/map_noFL"""
        return op.join(self.out_dir, "map_noFL")

    @property
    def quivered_dir(self):
        """Return $root_dir/quivered"""
        return op.join(self.root_dir, "quivered")

    @property
    def quivered_log_dir(self):
        """Return $log_dir/quivered"""
        return op.join(self.log_dir, "quivered")

    @property
    def nfl_all_pickle_fn(self):
        """Return $root_dir/$nfl_dir/nfl.all.partial_uc.pickle,
        this pickle file has all the paitial uc."""
        return op.join(self.nfl_dir, "nfl.all.partial_uc.pickle")

    @property
    def final_consensus_fa(self):
        """Return final consensus Fasta file."""
        return op.join(self.out_dir, "final.consensus.fasta")

    @property
    def final_consensus_sa(self):
        """Return suffix array of the final consensus Fa file."""
        return op.join(self.out_dir, "final.consensus.fasta.sa")

    @property
    def final_dazz_db(self):
        """Return final.consensus.dazz.fasta.db"""
        return op.join(self.out_dir, "final.consensus.dazz.fasta.db")

    @property
    def final_pickle_fn(self):
        """Return $root_dir/output/final.pickle"""
        return op.join(self.out_dir, "final.pickle")

    @property
    def submitted_quiver_jobs_log(self):
        """Return $root_dir/log/submitted_quiver_jobs.txt"""
        return op.join(self.log_dir, 'submitted_quiver_jobs.txt')

    def nfl_fa_i(self, i):
        """Return the i-th splitted chunk of nfl reads.
           $root_dir/output/map_noFL/input.split_{0:03d}.fa

           NOTE: make sure this agrees with io.FastaSplitter._out_fn()
        """
        fa_name = "input.split_{0:03d}.fasta".format(i)
        return op.join(self.nfl_dir, fa_name)

    def nfl_pickle_i(self, i):
        """Return the picke file of the i-th chunk of nfl reads.
           $root_dir/output/map_noFL/input.split_{0:03d}.fa.partial_uc.pickle
        """
        return self.nfl_fa_i(i) + ".partial_uc.pickle"

    def nfl_done_i(self, i):
        """Return the done file of the i-th chunk of nfl reads.
           $root_dir/output/map_noFL/input.split_{0:03d}.fa.partial_uc.pickle.DONE
        """
        return self.nfl_pickle_i(i) + ".DONE"

    def nfl_script_i(self, i):
        """Return the 'ice_partial' script file of the i-th chunk of nfl reads.
           $root_dir/scripts/input.split_{0:03d}.fa.partial_uc.sh
        """
        return op.join(self.script_dir,
                       op.basename(self.nfl_pickle_i(i)) + ".partial_uc.sh")

    def cluster_dir(self, cid):
        """Return directory path for the i-th cluster, i in [0,...]"""
        return op.join(self.tmp_dir, str(int(cid) / 10000), 'c' + str(cid))

    def raw_fa_of_cluster(self, cid):
        """Return $cluster_dir/in.raw_with_partial.fa, which
        contains the unrolled sequence of zmws belonging to cluster
        cid."""
        return op.join(self.cluster_dir(cid), "in.raw_with_partial.fasta")

    def g_consensus_fa_of_cluster(self, cid):
        """Return $cluster_dir(cid)/g_consensus.fa.
        Whenever this is changed, the ice_pbdagcon.py command needs
        to be changed accordingly.
        """
        return op.join(self.cluster_dir(cid), "g_consensus.fasta")

    def g_consensus_ref_fa_of_cluster(self, cid):
        """Return $cluster_dir(cid)/g_consensus_ref.fa.
        Whenever this is changed, the ice_pbdagcon.py command needs
        to be changed accordingly.
        """
        return op.join(self.cluster_dir(cid), "g_consensus_ref.fasta")

    def first_seq_fa_of_cluster(self, cid):
        """Return $cluster_dir(cid)/in.fa.1stseq.fa"""
        return op.join(self.cluster_dir(cid), "in.fa.1stseq.fasta")

    def sam_of_cluster(self, cid):
        """Return $cluster_dir/out.sam.
        This is a a sam file produced by blasr, aligning unrolled sequences
        of reads belonging to this cluster (e.g., in.raw_with_partial.fa)
        to consensus sequence of this cluster
        (e.g. g_consensus.fa or
              g_consensus_ref.fa or
              in.fa.1stseq.fa).
        """
        return op.join(self.cluster_dir(cid), "out.sam")

    def add_log(self, msg, level=logging.DEBUG):
        """Add a message to log_f and logging.info|debug|error."""
        msg = "[" + str(self.prog_name) + "] " + msg
        level_name = ""
        if level == logging.INFO:
            level_name = "[INFO]"
            logging.info(msg)
        elif level == logging.ERROR:
            level_name = "[ERROR]"
            logging.error(msg)
        elif level == logging.WARNING:
            level_name = "[WARN]"
            logging.warn(msg)
        else:
            level_name = "[DEBUG]"
            logging.debug(msg)

        if not self.no_log_f:
            self.log_f.write("{ln} [{pn}] {t}: {msg}\n".format(
                ln=level_name, pn=self.prog_name, t=now_str(), msg=msg))

    def close_log(self):
        """Close log file before exit."""
        if not self.no_log_f:
            self.add_log(self.prog_name + " completed.")
            self.log_f.close()

    def qsub_cmd_and_log(self, cmd):
        """Qsub the given command and write to log, raise a RunTimeError if
        failed to qsub, return qsub job id.

        The error message to display should look like:
            Failed to qsub CMD to SGE: {cmd}, {msg}\n
        """
        #msg = "Submitting CMD: {cmd}".format(cmd=cmd)
        #self.add_log(msg)
        _out, _code, _msg = backticks(cmd)
        if _code != 0:
            errMsg = "Failed to qsub CMD: {cmd}, {msg}.".\
                format(cmd=cmd, msg=_msg)
            self.add_log(errMsg, level=logging.ERROR)
            raise RuntimeError(errMsg)
        # Your job 596028 ("a.sh") has been submitted
        return str(_out).split()[2]

    def run_cmd_and_log(self, cmd, olog="", elog="", description=""):
        """Run the given command locally and write to log, raise a
        RunTimeError if failed to finish the job.
        olog: output log
        elog: error log

        The error message to display should look like:
            CMD exited with a non-zero code: {cmd}, {msg}\n
            {description}\n
            Error log: {elog}\n
        """
        #msg = "Running CMD: {cmd}".format(cmd=cmd)
        #self.add_log(msg)
        _out, _code, _msg = backticks(cmd)
        if _code != 0:
            errMsgs = ["CMD exited with a non-zero code: {cmd}, {msg}".
                       format(cmd=cmd, msg=_msg)]
            if len(description) != 0:
                errMsgs.append("{description}".format(description=description))
            if len(elog) != 0:
                errMsgs.append("Error log: {elog}".format(elog=elog))
            if len(olog) != 0:
                errMsgs.append("Out log: {olog}".format(olog=olog))
            errMsg = "\n".join(errMsgs)
            self.add_log(errMsg, level=logging.ERROR)
            raise RuntimeError(errMsg)

    def write_report(self, report_fn, uc, partial_uc=None):
        """
        Write a CSV report to report_fn, each line contains three columns:
            cluster_id, read_id and read_type
        """
        self.add_log("Writing a csv report of cluster -> FL{nfl} reads to {f}".
                     format(nfl="/NonFL" if partial_uc is not None else "",
                            f=report_fn), level=logging.INFO)
        with open(report_fn, 'w') as f:
            f.write("cluster_id,read_id,read_type\n")
            for c in uc.keys():
                for r in uc[c]:
                    f.write("c{c},{r},FL\n".format(r=r, c=c))
                if partial_uc is not None:
                    for r in partial_uc[c]:
                        f.write("c{c},{r},NonFL\n".format(r=r, c=c))

    def write_summary(self, summary_fn, isoforms_fa, hq_fa=None, lq_fa=None):
        """Extract number of consensus isoforms predicted, and total
        number of bases in all consensuus isoforms from isoforms_fa and write
        the two attributes to summary_fn.

        if hq_fa (polished high-quality isoforms) is not None, report
            the number of polished hq clusters
        if lq_fa (polished high-quality isoforms) is not None, report
            the number of polished hq clusters
        """
        self.add_log("Writing a summary to {f}".format(f=summary_fn),
                     level=logging.INFO)
        try:
            summary = ClusterSummary()

            with FastaReader(isoforms_fa) as reader:
                for r in reader:
                    summary.numConsensusIsoforms += 1
                    summary.numTotalBases += len(r.sequence)

            if hq_fa is not None:
                summary.num_polished_hq_isoforms = 0
                with FastaReader(hq_fa) as reader:
                    for r in reader:
                        summary.num_polished_hq_isoforms += 1
            if lq_fa is not None:
                summary.num_polished_lq_isoforms = 0
                with FastaReader(lq_fa) as reader:
                    for r in reader:
                        summary.num_polished_lq_isoforms += 1
            summary.write(summary_fn)
        except ZeroDivisionError:
            errMsg = "No consensus isoforms predicted."
            self.add_log(errMsg, level=logging.ERROR)
            raise RuntimeError(errMsg)


def wait_for_sge_jobs_worker(cmd):
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        errMsg = "Failed to qsub CMD: {cmd}, {msg}.".\
            format(cmd=cmd, msg=_msg)
        raise RuntimeError(errMsg)
    # Your job 596028 ("a.sh") has been submitted
    return str(_out).split()[2]

def wait_for_sge_jobs(cmd, jids, timeout):
    """
    This replaces the original qsub -sync y -hold_jid j1,j2..... command
    which can still be hung if certain jobs got stuck.

    If timeout occurs, simply qdel all jids (ignoring whether they exist or not)
    and let the main function that calls it handle what to do
    """
    def get_active_jids():
        stuff = os.popen("qstat").read().strip().split('\n')
        for x in stuff[2:]:
            job_id = x.split()[0]
            yield job_id

    p = Process(target=wait_for_sge_jobs_worker, args=(cmd,))
    p.start()
    p.join(timeout)
    if p.is_alive(): # timed out
        active_jids = [x for x in get_active_jids()]
        while len(active_jids) > 0:
            for jid in active_jids:
                kill_cmd = "qdel " + str(jid)
                backticks(kill_cmd) # don't care whether it worked
            time.sleep(3) # wait 3 sec for qdel to take effect....
            active_jids = [x for x in get_active_jids()] # make sure qdel really worked
        return "TIMEOUT"
    return "SUCCESS"


