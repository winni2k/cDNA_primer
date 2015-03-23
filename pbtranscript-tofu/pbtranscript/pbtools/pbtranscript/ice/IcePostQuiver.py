#!/usr/bin/env python

"""Post quiver, pick up the best conesnsus isoforms."""
import os, sys, re
import logging
import os.path as op
from collections import defaultdict
from cPickle import load
from time import sleep
from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.PBTranscriptOptions import \
        add_polished_isoforms_arguments
from pbtools.pbtranscript.Utils import phred_to_qv, \
        get_all_files_in_dir, ln
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbcore.io import FastaWriter, FastqReader, FastqWriter


class IcePostQuiver(IceFiles):
    """check if quiver jobs are finished and quiver results are compeleted.
       If quiver jobs are completed, pick up high QV consensus isoforms.
       If use_sge is True and quiver jobs are still running,
           * If quit_if_not_done is True, exit.
           * If quit_if_not_done is False, wait till quiver jobs are finished.

    """
    # ignore first 200 bp on 5' end...
    #   Quiver usually can't  call it very well + less coverage
    qv_trim_5 = 200
    qv_trim_3 = 50 # """ignore last 50 bp on the 3' end"""
    # total number of allowed expected base errors within
    # seq[qv_trim_5:seq_trim_3]
    qv_max_err = 10

    def __init__(self, root_dir, hq_isoforms_fa=None, hq_isoforms_fq=None,
                 lq_isoforms_fa=None, lq_isoforms_fq=None,
                 use_sge=False, quit_if_not_done=True):
        super(IcePostQuiver, self).__init__(prog_name="IcePostQuiver",
              root_dir=root_dir)
        self.hq_isoforms_fa = hq_isoforms_fa
        self.hq_isoforms_fq = hq_isoforms_fq
        self.lq_isoforms_fa = lq_isoforms_fa
        self.lq_isoforms_fq = lq_isoforms_fq
        self.use_sge = use_sge
        self.quit_if_not_done = quit_if_not_done
        self.fq_filenames = []
        self._validate_inputs()

    def get_existing_binned_quivered_fq(self):
        """Return all existing quivered fq files for binned clusters."""
        pattern = r"c(\d+)to(\d+)"  # e.g. c0to214
        fs = get_all_files_in_dir(self.quivered_dir,
                                  extension="quivered.fq")
        return [f for f in fs if re.search(pattern, f) is not None]

    def _validate_inputs(self):
        """Validate if logs  and pickle for non-full-length reads exist."""
        errMsg = ""
        if not op.exists(self.nfl_all_pickle_fn):
            errMsg = "The pickle of all non-full-length reads " + \
                     "{f} does not exists.".format(f=self.nfl_all_pickle_fn)
        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise IOError(errMsg)

    def check_quiver_jobs_completion(self):
        """Check whether quiver jobs are completed.
        submitted_quiver_jobs.txt should have format like:
        <job_id> \t ./quivered/<range>.sh

        (1) if all jobs are done and files are there return True
        (2) if all jobs are done but some files incomplete ask if to resubmit
        (3) if not all jobs are done, just quit
        fq_filenames contains all the finished fastq files.
        """
        self.add_log("Checking if quiver jobs are completed.")
        done_flag = True
        bad_sh = []
        self.fq_filenames = []
        submitted = {}
        self.add_log("Submitted quiver jobs are at {f}:".
                     format(f=self.submitted_quiver_jobs_log))

        with open(self.submitted_quiver_jobs_log, 'r') as f:
            for line in f:
                a, b = line.strip().split('\t')
                if a == 'local':
                    submitted[b] = b
                else:
                    submitted[a] = b

        stuff = os.popen("qstat").read().strip().split('\n')
        # first two lines are header
        running_jids = []
        for x in stuff[2:]:
            job_id = x.split()[0]
            running_jids.append(job_id)
            if job_id in submitted:
                self.add_log("job {0} is still running.".format(job_id))
                done_flag = False

        for job_id, sh_name in submitted.iteritems():
            fq_filename = op.join(self.quivered_dir,
                op.basename(sh_name).replace('.sh', '.quivered.fq'))

            if not op.exists(fq_filename) or \
                os.stat(fq_filename).st_size == 0:
                if job_id in running_jids:  #still running, pass
                    done_flag = False
                else:
                    self.add_log("job {0} is completed but {1} is still empty!".
                                 format(job_id, fq_filename))
                    bad_sh.append(submitted[job_id])
            else:
                self.add_log("job {0} is done".format(job_id))
                self.fq_filenames.append(fq_filename)

        if not done_flag:
            if len(bad_sh) == 0:
                return "RUNNING"
            else:
                self.add_log("The following jobs were completed but " +
                    "no output file. Please check and resubmit: " +
                    "\n{0}\n".format('\n'.join(bad_sh)))
                return "FAILED"
        else:
            return "DONE"

    @property
    def quivered_good_fa(self):
        """Return $root_dir/all_quivered.hq.a_b_c.fasta"""
        return op.join(self.root_dir,
                       "all_quivered_hq.{a}_{b}_{c}.fasta".format(
                       a=self.qv_trim_5,
                       b=self.qv_trim_3,
                       c=self.qv_max_err))

    @property
    def quivered_good_fq(self):
        """Return $root_dir/all_quivered_hq.a_b_c.fq"""
        return op.join(self.root_dir,
                       "all_quivered_hq.{a}_{b}_{c}.fastq".format(
                       a=self.qv_trim_5,
                       b=self.qv_trim_3,
                       c=self.qv_max_err))

    @property
    def quivered_bad_fa(self):
        """Return $root_dir/all_quivered_lq.fa"""
        return op.join(self.root_dir, "all_quivered_lq.fasta")

    @property
    def quivered_bad_fq(self):
        """Return $root_dir/all_quivered_lq.fq"""
        return op.join(self.root_dir, "all_quivered_lq.fastq")

    def pickup_best_clusters(self, fq_filenames):
        """Pick up hiqh QV clusters."""
        self.add_log("Picking up the best clusters according to QVs from {fs}.".
                     format(fs=", ".join(fq_filenames)))
        a = load(open(self.final_pickle_fn))
        uc = a['uc']
        quivered = {}

        for fq in fq_filenames:
            self.add_log("Looking at quivered fq {f}".format(f=fq))
            for r in FastqReader(fq):
                cid = r.name.split('|')[0]
                if cid.endswith('_ref'):
                    cid = cid[:-4]
                cid = int(cid[1:])
                quivered[cid] = r

        good = []

        for cid, r in quivered.iteritems():
            q = [phred_to_qv(x) for x in r.quality]
            if sum(q[self.qv_trim_5: -self.qv_trim_3]) <= self.qv_max_err:
                good.append(cid)

        partial_uc = load(open(self.nfl_all_pickle_fn))['partial_uc']
        partial_uc2 = defaultdict(lambda: [])
        partial_uc2.update(partial_uc)

        self.add_log("Writing hiqh-quality isoforms to {f}|fq".
                     format(f=self.quivered_good_fa))
        self.add_log("Writing low-quality isoforms to {f}|fq".
                     format(f=self.quivered_bad_fa))
        with FastaWriter(self.quivered_good_fa) as good_fa_writer, \
             FastaWriter(self.quivered_bad_fa) as bad_fa_writer, \
             FastqWriter(self.quivered_good_fq) as good_fq_writer, \
             FastqWriter(self.quivered_bad_fq) as bad_fq_writer:
            for cid in quivered:
                r = quivered[cid]
                newname = "c{cid}/f{flnc_num}p{nfl_num}/{read_len}".\
                        format(cid=cid,
                               flnc_num=len(uc[cid]),
                               nfl_num=len(partial_uc2[cid]),
                               read_len=len(r.sequence))

                if cid in good:
                    self.add_log("processing quivered cluster {c} --> good.".
                                 format(c=cid))
                    good_fa_writer.writeRecord(newname, r.sequence)
                    good_fq_writer.writeRecord(newname, r.sequence, r.quality)
                else:
                    self.add_log("processing quivered cluster {c} --> bad.".
                                 format(c=cid))
                    bad_fa_writer.writeRecord(newname, r.sequence)
                    bad_fq_writer.writeRecord(newname, r.sequence, r.quality)

        self.add_log("-" * 60, level=logging.INFO)
        self.add_log("High-quality Quivered consensus written " +
                     "to:\n{0}\n{1}\n".format(self.quivered_good_fa,
                     self.quivered_good_fq))
        self.add_log("Low-qulality Quivered consensus written " +
                     "to:\n{0}\n{1}".format(self.quivered_bad_fa,
                     self.quivered_bad_fq))
        self.add_log("-" * 60, level=logging.INFO)

    def run(self):
        """Check all quiver jobs are running, failed or done. Write high-quality
        consensus and low-quality consensus to all_quivered.good|bad.fa|fq.
        """
        job_stats = self.check_quiver_jobs_completion()
        self.add_log("quiver job status: {s}".format(s=job_stats))

        if self.use_sge is not True and job_stats != "DONE":
            self.add_log("quiver jobs were not submitted via sge, " +
                         "however are still incomplete. Please check.",
                         level=logging.ERROR)
            return -1
        elif self.use_sge is True:
            while job_stats != "DONE":
                self.add_log("Sleeping for 180 seconds.")
                sleep(180)
                job_stats = self.check_quiver_jobs_completion()
                if job_stats == "DONE":
                    break
                elif job_stats == "FAILED":
                    self.add_log("There are some failed jobs. Please check.",
                                 level=logging.ERROR)
                    return 1
                elif job_stats == "RUNNING":
                    self.add_log("There are jobs still running, waiting...",
                                 level=logging.INFO)
                    if self.quit_if_not_done is True:
                        return 0
                else:
                    msg = "Unable to recognize job_stats {s}".format(job_stats)
                    self.add_log(msg, logging.ERROR)
                    raise ValueError(msg)

        self.pickup_best_clusters(self.fq_filenames)

        self.add_log("Creating polished high quality consensus isoforms.")
        if self.hq_isoforms_fa is not None:
            ln(self.quivered_good_fa, self.hq_isoforms_fa)
        if self.hq_isoforms_fq is not None:
            ln(self.quivered_good_fq, self.hq_isoforms_fq)

        self.add_log("Creating polished low quality consensus isoforms.")
        if self.lq_isoforms_fa is not None:
            ln(self.quivered_bad_fa, self.lq_isoforms_fa)
        if self.lq_isoforms_fq is not None:
            ln(self.quivered_bad_fq, self.lq_isoforms_fq)

        self.close_log()


def set_parser(parser):
    """Set argument parser."""
    parser.add_argument("root_dir",
                        type=str,
                        help="Root directory for ICE, " + \
                             "e.g., SMRTPortal data/clusterOutDir")
    parser = add_polished_isoforms_arguments(parser)
    parser.add_argument("--use_sge",
                        default=False,
                        dest="use_sge",
                        action="store_true",
                        help="quiver jobs have been submitted to sge."
                             "Check qstat")
    parser.add_argument("--quit_if_not_done",
                        default=False,
                        dest="quit_if_not_done",
                        action="store_true",
                        help="Quit if quiver jobs haven't been completed.")


from pbcore.util.ToolRunner import PBToolRunner
class IcePostQuiverRunner(PBToolRunner):
    """IcePostQuiver runner"""
    def __init__(self):
        desc = "Post-quiver processing, selecting high QV consensus isoforms."
        PBToolRunner.__init__(self, desc)
        set_parser(self.parser)

    def getVersion(self):
        """Get version string"""
        return get_version()

    def run(self):
        """Run"""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))
        args = self.args
        try:
            obj = IcePostQuiver(root_dir=args.root_dir,
                                hq_isoforms_fa=args.hq_isoforms_fa,
                                hq_isoforms_fq=args.hq_isoforms_fq,
                                lq_isoforms_fa=args.lq_isoforms_fa,
                                lq_isoforms_fq=args.lq_isoforms_fq,
                                use_sge=args.use_sge,
                                quit_if_not_done=args.quit_if_not_done)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = IcePostQuiverRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

