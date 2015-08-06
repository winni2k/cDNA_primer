#!/usr/bin/env python

"""Call quiver to polish consensus isoforms created by ICE."""

import os.path as op
import time
from datetime import datetime
import random
import logging
import shutil
from cPickle import load
from math import ceil
from collections import defaultdict
from multiprocessing.pool import ThreadPool
from pbtools.pbtranscript.PBTranscriptOptions import  add_fofn_arguments, \
    add_sge_arguments, add_cluster_root_dir_as_positional_argument
from pbtools.pbtranscript.Utils import mkdir, real_upath, nfs_exists
from pbtools.pbtranscript.ice.IceUtils import get_the_only_fasta_record, \
    get_files_from_fofn, is_blank_sam, concat_sam, \
    blasr_sam_for_quiver, write_in_raw_fasta, write_in_raw_fasta_starhelper
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.io.FastaRandomReader import MetaSubreadFastaReader


class IceQuiver(IceFiles):

    """Ice Quiver."""

    desc = "After assigning all non-full-length reads to unpolished " + \
           "consensus isoforms created by ICE, call quiver to polish " + \
           "these isoforms."

    def __init__(self, root_dir, bas_fofn, fasta_fofn, sge_opts,
                 prog_name=None):
        # Initialize super class IceFiles.
        prog_name = "IceQuiver" if prog_name is None else prog_name
        super(IceQuiver, self).__init__(prog_name=prog_name,
                                        root_dir=root_dir, bas_fofn=bas_fofn,
                                        fasta_fofn=fasta_fofn)
        self.sge_opts = sge_opts

    def validate_inputs(self):
        """Validate input fofns, and root_dir, log_dir, tmp_dir,
        create quivered_dir and quivered_log_dir"""
        self.add_log("Validating inputs.")

        # Create directories: root_dir/quivered and root_dir/log_dir/quivered
        try:
            mkdir(self.quivered_dir)
            mkdir(self.quivered_log_dir)
        except OSError:
            # Multiple ice_quiver_i jobs may run at the same time and try to
            # mkdir, race condition may happen, so ignore OSError here.
            pass

        errMsg = ""

        if not nfs_exists(self.log_dir) or not op.isdir(self.log_dir):
            errMsg = "Log dir {l} is not an existing directory.".\
                format(l=self.log_dir)
        elif self.bas_fofn is None:
            errMsg = "Please specify bas_fofn (e.g. input.fofn)."
        elif not nfs_exists(self.bas_fofn):
            errMsg = "bas_fofn {f} ".format(f=self.bas_fofn) + \
                     "which contains bas/bax.h5 files does not exist."
        elif self.fasta_fofn is None:
            errMsg = "Please make sure ice_make_fasta_fofn has " + \
                     "been called, and specify fasta_fofn."
        elif not nfs_exists(self.fasta_fofn):
            errMsg = "Input fasta_fofn {f} does not exists.".\
                     format(f=self.fasta_fofn)
            fasta_files = get_files_from_fofn(self.fasta_fofn)
            for fasta_file in fasta_files:
                if not nfs_exists(fasta_file):
                    errMsg = "A file {f} in fasta_fofn does not exist.".\
                             format(f=fasta_file)
        elif not nfs_exists(self.nfl_all_pickle_fn):
            #"output/map_noFL/noFL.ALL.partial_uc.pickle"):
            errMsg = "Pickle file {f} ".format(f=self.nfl_all_pickle_fn) + \
                     "which assigns all non-full-length reads to isoforms " + \
                     "does not exist. Please check 'ice_partial.py *' are " + \
                     "all done."
        elif not nfs_exists(self.final_pickle_fn):
            errMsg = "Pickle file {f} ".format(f=self.final_pickle_fn) + \
                     "which assigns full-length non-chimeric reads to " + \
                     "isoforms does not exist."

        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise IOError(errMsg)

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
        return self._quivered_bin_prefix(first, last) + ".ref.fa"

    def cmph5_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.cmp.h5"""
        return self._quivered_bin_prefix(first, last) + ".cmp.h5"

    def fq_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.quivered.fq
        this is quivered fq output. Whenever this is changed, change
        IceQuiverPostprocess accordingly.
        """
        return self._quivered_bin_prefix(first, last) + ".quivered.fq"

    def script_of_quivered_bin(self, first, last):
        """Return $_quivered_bin_prefix.sh"""
        return self._quivered_bin_prefix(first, last) + ".sh"

    def create_raw_fas_for_clusters_in_bin(self, cids, d, uc, partial_uc):
        """
        Create raw subreads fasta files for clusters in cids.
        For each cluster k in cids,
        * Collect raw subreads of zmws associated with cluster k
          in either uc or partial_uc.

        cids --- cluster ids
        d --- MetaSubreadsFastaReader
        uc --- uc[k] returns fl ccs reads associated with cluster k
        partial_uc --- partial_uc[k] returns nfl ccs reads associated with cluster k
        
        (Liz) for Quiver, subsample down to max 100 (with uc having priority over partial_uc)
        """
        #data_queue = []
        #fake_func = lambda x: x
        for k in cids:  # for each cluster k

            # $root_dir/tmp/?/c{k}/in.raw_with_partial.fa
            raw_fa = self.raw_fa_of_cluster(k)

            in_seqids = uc[k]
            if len(in_seqids) > 100:
                in_seqids = random.sample(in_seqids, 100)
            else:
                in_seqids += random.sample(partial_uc[k], min(len(partial_uc[k]), 100 - len(in_seqids)))
            # write cluster k's associated raw subreads to raw_fa
            #data_queue.append([d, in_seqids, raw_fa, True])
            write_in_raw_fasta(input_fasta_d=d,
                               in_seqids=in_seqids,
                               out_fa=raw_fa,
                               ignore_keyerror=True)
        #p = ThreadPool(processes=3)
        #rets = p.map(write_in_raw_fasta_starhelper, data_queue)
        #p.close()
        #p.join()

    def create_sams_for_clusters_in_bin(self, cids, refs):
        """
        Create sam files for clusters in cids.
        For each cluster k in cids,
        * Call blasr to align its associated subreads in raw_fa to its consensus
          sequence as reference.

        cids --- cluster ids
        refs --- refs[k] -> consensus seq of cluster k

        This function has to be called after prepare_raw_fa_for_clusters is done.

        """
        for k in cids:  # for each cluster k

            # $root_dir/tmp/?/c{k}/in.raw_with_partial.fa
            raw_fa = self.raw_fa_of_cluster(k)

            if not op.exists(raw_fa):
                raise IOError("raw_fa {f} does not exist. ".format(f=raw_fa) +
                              "Please check raw fasta of this bin is created.")

            blasr_sam_for_quiver(
                input_fasta=raw_fa,
                ref_fasta=refs[k],
                out_sam_filename=self.sam_of_cluster(k),
                run_cmd=True,
                blasr_nproc=self.sge_opts.blasr_nproc)

    def concat_valid_sams_and_refs_for_bin(self, cids, refs):
        """
        Concat sam files and reference sequences of all valid clusters
        in bin to create a big sam and a big ref.
        A cluser is not valid if (1) or (2)
            (1) identical sequences already exists in another cluster
                (rare, but happens)
            (2) the alignment is empty (also rare, but happens)
        Return valid_cids, a list of valid cluster ids
        """
        first, last = cids[0], cids[-1]
        bin_sam_file = self.sam_of_quivered_bin(first, last)
        bin_ref_fa = self.ref_fa_of_quivered_bin(first, last)

        self.add_log("Concatenating reference files between " +
                     "{first} and {last}.".format(first=first, last=last))
        valid_sam_files = []
        valid_cids = []
        seqs_seen = {}
        with open(bin_ref_fa, 'w') as bin_ref_fa_writer:
            for cid in cids:
                fname = self.sam_of_cluster(cid)
                if not is_blank_sam(fname):
                    ref_rec = get_the_only_fasta_record(refs[cid])
                    name = ref_rec.name.strip()
                    #if '/' in name:
                    #    # convert both c{cid} and c{cid}/0_len to c{cid}
                    #    name = name[:name.find('/')]
                    seq = ref_rec.sequence.strip()
                    if seq not in seqs_seen:
                        valid_sam_files.append(fname)
                        valid_cids.append(cid)
                        seqs_seen[seq] = cid
                        # concate valid ref files, avoid 'cat ...' hundreds
                        # or even thousands of files due to linux cmd line
                        # length limits
                        bin_ref_fa_writer.write(">{0}\n{1}\n".
                                                format(name, seq))
                    else:
                        self.add_log("ignoring {0} because identical sequence!".format(cid))
                else:
                    self.add_log("ignoring {0} because no alignments!".format(cid))

        if len(valid_sam_files) == 0:
            self.add_log("No alignments were found for clusters between " +
                         "{first} and {last}.".format(first=first, last=last),
                         level=logging.WARNING)
            assert(len(valid_cids) == 0)
        else:
            self.add_log("Concatenating sam files between " +
                         "{first} and {last}.".format(first=first, last=last))
            # concat valid sam files
            concat_sam(valid_sam_files, bin_sam_file)
            self.add_log("Concatenation done")

        return valid_cids

    def quiver_cmds_for_bin(self, cids, quiver_nproc=2):
        """Return a list of quiver related cmds, to convert sam & ref to cmp.h5
        and call quiver, including samtoh5, loadPulses, comph5tools.py,
        samtools, loadChemistry, quiver...
        """
        first, last = cids[0], cids[-1]
        self.add_log("Creating quiver cmds for c{first} to c{last}".
                     format(first=first, last=last))

        bin_sam_file = self.sam_of_quivered_bin(first, last)
        bin_ref_fa = self.ref_fa_of_quivered_bin(first, last)
        bin_cmph5 = self.cmph5_of_quivered_bin(first, last)
        bin_fq = self.fq_of_quivered_bin(first, last)

        cmds = []
        cmds.append("samtoh5 {sam} {ref} {cmph5} -smrtTitle".format(
            sam=real_upath(bin_sam_file),
            ref=real_upath(bin_ref_fa),
            cmph5=real_upath(bin_cmph5)))
        # (Liz) don't gzip the sa
        #cmds.append("gzip {sam}".format(sam=real_upath(bin_sam_file)))
        metrics = ["QualityValue", "InsertionQV", "MergeQV", "DeletionQV",
                   "DeletionTag", "SubstitutionTag", "SubstitutionQV"]
        cmds.append("loadPulses {bas_fofn} ".
                    format(bas_fofn=real_upath(self.bas_fofn)) +
                    "{cmph5} ".format(cmph5=real_upath(bin_cmph5)) +
                    "-byread -metrics " + ",".join(metrics))
        cmds.append("cmph5tools.py sort {cmph5}".
                    format(cmph5=real_upath(bin_cmph5)))
        cmds.append("samtools faidx {ref}".format(ref=real_upath(bin_ref_fa)))
        cmds.append("loadChemistry.py {bas_fofn} {cmph5}".
                    format(bas_fofn=real_upath(self.bas_fofn),
                           cmph5=real_upath(bin_cmph5)))
        cmds.append("quiver {cmph5} ".format(cmph5=real_upath(bin_cmph5)) +
                    "-v -j{n} ".format(n=quiver_nproc) +
                    "-r {ref} ".format(ref=real_upath(bin_ref_fa)) +
                    "-o {fq}".format(fq=real_upath(bin_fq)))
        return cmds

    def create_quiver_sh_for_bin(self, cids, cmds):
        """
        Write quiver cmds to a bash script, e.g., quivered/c{}to{}.sh,
        return script file path.
        """
        first, last = cids[0], cids[-1]
        bin_sh = self.script_of_quivered_bin(first, last)
        self.add_log("Creating quiver bash script {f} for c{first} to c{last}.".
                     format(f=bin_sh, first=first, last=last))
        with open(bin_sh, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("\n".join(cmds))

        return bin_sh

    def submit_todo_quiver_jobs(self, todo, submitted, sge_opts):
        """
        todo --- a list of sh scripts to run
        submitted --- a list of sh scripts which have been submitted
        sge_opts --- SGE options, including
                     use_sge, whether or not to use sge
                     max_sge_jobs, maximum number sge jobs to submit
                     quiver_nproc, number of nproc per job
                     unique_id, unique id to name qsub jobs
        """
        time0 = datetime.now()
        if sge_opts.use_sge is not True or \
           sge_opts.max_sge_jobs == 0:  # don't use SGE
            for job in todo:
                elog = op.join(self.quivered_log_dir,
                               op.basename(job) + ".elog")
                olog = op.join(self.quivered_log_dir,
                               op.basename(job) + ".olog")
                cmd = "bash " + real_upath(job) + " 1>{olog} 2>{elog}".\
                      format(olog=real_upath(olog), elog=real_upath(elog))
                self.run_cmd_and_log(cmd, olog=olog, elog=elog,
                                     description="Failed to run Quiver")
                submitted.append(("local", job))
            todo = []
        else:
            while len(todo) > 0:
                n = min(sge_opts.max_sge_jobs, len(todo))
                for job in todo[:n]:
                    # ex: Your job 8613116 ("c20to70.sh") has been submitted
                    elog = op.join(self.quivered_log_dir,
                                   op.basename(job) + ".elog")
                    olog = op.join(self.quivered_log_dir,
                                   op.basename(job) + ".olog")
                    jid = "ice_quiver_{unique_id}_{name}".format(
                        unique_id=self.sge_opts.unique_id,
                        name=op.basename(job))

                    qsub_cmd = "qsub"
                    if self.sge_opts.sge_queue is not None:
                        qsub_cmd += " -q " + self.sge_opts.sge_queue
                    qsub_cmd += " -pe {env} {n} ".format(n=sge_opts.quiver_nproc, env=sge_opts.sge_env_name) + \
                               "-cwd -S /bin/bash -V " + \
                               "-e {elog} ".format(elog=real_upath(elog)) +\
                               "-o {olog} ".format(olog=real_upath(olog)) +\
                               "-N {jid} ".format(jid=jid) + \
                               "{job}".format(job=real_upath(job))
                    job_id = self.qsub_cmd_and_log(qsub_cmd)

                    submitted.append((job_id, job))
                    todo.remove(job)
                # end of for job in todo[:n]
            # end of while len(todo) > 0
        # end of else (use sge)
        self.add_log("Total time submitting todo quiver jobs: {0}".format(datetime.now()-time0))

    def create_a_quiver_bin(self, cids, d, uc, partial_uc, refs, sge_opts):
        """Put clusters in cids together into a bin. In order to polish
        consensus of clusters in the bin, prepare inputs and create a quiver
        bash script to run later.

        (1) For each cluster k in cids, obtain unrolled sequences of all zmws
            belonging to this cluster, and save in raw_fa_of_cluster(k)
        (2) For each cluster k in cids, call blasr to align raw_fa_of_cluster to
            its consensus sequence and create sam_of_cluster(k).
        (3) Concat all sam files of `valid` clusters to sam_of_quivered_bin, and
            concat ref seqs of all `valid` clusters to ref_fa_of_quivered_bin
        (4) Make commands including
                samtoh5, loadPulses, cmph5tools.py, loadChemistry, ..., quiver
            in order to convert sam_of_quivered_bin to cmph5_of_quivered_bin.
            Write these commands to script_of_quivered_bin
              * qsub all jobs later when scripts of all quivered bins are done.
              * or execute scripts sequentially on local machine
        """
        # For each cluster in bin, create its raw subreads fasta file.
        time0 = datetime.now()
        self.create_raw_fas_for_clusters_in_bin(cids=cids, d=d, uc=uc,
                                                partial_uc=partial_uc)
        self.add_log("Total time for create_raw_fas_for_clusters_in_bin: {0}".format(datetime.now()-time0))

        # For each cluster in bin, align its raw subreads to ref to build a sam
        self.create_sams_for_clusters_in_bin(cids=cids, refs=refs)

        # Concatenate sam | ref files of 'valid' clusters in this bin to create
        # a big sam | ref file.
        valid_cids = self.concat_valid_sams_and_refs_for_bin(cids=cids, refs=refs)

        # quiver cmds for this bin
        cmds = []
        if len(valid_cids) != 0:
            cmds = self.quiver_cmds_for_bin(cids=cids,
                                            quiver_nproc=sge_opts.quiver_nproc)
        else:
            cmds = ["echo no valid clusters in this bin, skip..."]

        # Write quiver cmds for this bin to $root_dir/quivered/c{}_{}.sh
        return self.create_quiver_sh_for_bin(cids=cids, cmds=cmds)

    def create_quiver_bins(self, d, uc, partial_uc, refs, keys, start, end,
                           sge_opts):
        """
        Create quiver bins by putting every 100 clusters into a bin.
        For each bin, create a bash script (e.g., script_of_quivered_bin).
        Return a list of scripts to run.
        """
        bin_scripts = []
        for i in xrange(start, end, 100):  # Put every 100 clusters to a bin
            cids = keys[i:min(end, i + 100)]
            bin_sh = self.create_a_quiver_bin(cids=cids, d=d, uc=uc,
                                              partial_uc=partial_uc,
                                              refs=refs, sge_opts=sge_opts)
            bin_scripts += bin_sh
        return bin_scripts

    def create_quiver_bins_and_submit_jobs(self, d, uc, partial_uc, refs, keys,
                                           start, end, submitted, sge_opts):
        """
        Put every 100 clusters together and create bins. Create a bash script
        (e.g., script_of_quivered_bin), for each bin, and submit the script
        either using qsub or running it locally.
        return all bash scripts in a list.
        """
        all_todo = []
        for i in xrange(start, end, 100):  # Put every 100 clusters to a bin
            cids = keys[i:min(end, i + 100)]
            time0 = time.time()
            bin_sh = self.create_a_quiver_bin(cids=cids, d=d, uc=uc,
                                              partial_uc=partial_uc,
                                              refs=refs, sge_opts=sge_opts)
            all_todo += bin_sh
            self.add_log("DEBUG: Total time for create_a_quiver_bin: {0}".format(time.time()-time0))
            # assert bin_sh == self.script_of_quivered_bin(first, last)
            # submit the created script of this quiver bin
            time1 = time.time()
            self.submit_todo_quiver_jobs(todo=[bin_sh], submitted=submitted,
                                         sge_opts=sge_opts)
            self.add_log("DEBUG: Total time for submit_todo_quiver_jobs: {0}".format(time.time()-time1))
        # end of for i in xrange(start, end, 100):
        return all_todo

    @property
    def report_fn(self):
        """Return a csv report with cluster_id, read_id, read_type."""
        return op.join(self.out_dir, "cluster_report.FL_nonFL.csv")

    def load_pickles(self):
        """Load uc and refs from final_pickle_fn, load partial uc from
        nfl_all_pickle_fn, return (uc, partial_uc. refs).
        """
        self.add_log("Loading uc from {f}.".format(f=self.final_pickle_fn))
        a = load(open(self.final_pickle_fn))
        uc = a['uc']
        refs = a['refs']

        self.add_log("Loading partial uc from {f}.".
                     format(f=self.nfl_all_pickle_fn))
        partial_uc = load(open(self.nfl_all_pickle_fn))['partial_uc']
        partial_uc2 = defaultdict(lambda: [])
        partial_uc2.update(partial_uc)
        return (uc, partial_uc2, refs)

    def index_fasta(self):
        """index subreads in fasta_fofn, return"""
        files = get_files_from_fofn(self.fasta_fofn)
        msg = "Indexing {0} fasta files, please wait.".format(len(files))
        self.add_log(msg)

        d = MetaSubreadFastaReader(files)
        self.add_log("Fasta files indexing done.")
        return d

    def submitted_quiver_jobs_log_of_chunk_i(self, i, num_chunks):
        """A txt file to save all submitted quiver jobs of the
        (i / num_chunks)-th workload. Format:
            job_id\tscript_path
        Return $root_dir/log/submitted_quiver_jobs.{i}of{num_chunks}.txt
        """
        return op.join(self.log_dir, "submitted_quiver_jobs.{i}of{N}.txt".
                                     format(i=i, N=num_chunks))

    def quiver_jobs_sh_of_chunk_i(self, i, num_chunks):
        """A bash file to save all quiver jobs of the
        (i / num_chunks)-th workload.
        Return $root_dir/log/quiver_jobs.{i}of{N}.sh
        """
        return op.join(self.log_dir, "quiver_jobs.{i}of{N}.sh".
                                     format(i=i, N=num_chunks))

    def process_chunk_i(self, i, num_chunks):
        """
        In order to distribute IceQuiver jobs by SMRTPipe using a fixed
        number of nodes, we divide quiver jobs into num_chunks workloads
        of roughly the same size, and are processing the i-th workload
        now.
        (1) load uc, partial_uc and refs from pickles and index subreads
            in fasta and save to d
        (2) write report if this is the first chunk (e.g, i==0)
        (3) Assume clusters are divided into num_chunks parts, process
            the i-th part.
        (4)

        """
        if (i >= num_chunks):
            raise ValueError("Chunk index {i} should be less than {N}.".
                             format(i=i, N=num_chunks))

        # load uc, partial_uc and refs from pickles,
        uc, partial_uc, refs = self.load_pickles()

        # Write report to quivered/cluster_report.FL_nonFL.csv
        if i == 0:
            self.write_report(report_fn=self.report_fn,
                              uc=uc, partial_uc=partial_uc)

        # index fasta files in fasta_fofn, and save to d
        d = self.index_fasta()

        # good = [x for x in uc if len(uc[x]) > 1 or len(partial_uc2[x]) >= 10]
        # bug 24984, call quiver on everything, no selection is needed.
        keys = sorted([x for x in uc])  # sort cluster ids

        # Compute number of clusters in i-th chunk
        num_clusters_per_chunk = int(ceil(len(keys) / float(num_chunks)))
        num_clusters_in_chunk_i = max(0, min(len(keys) - i * num_clusters_per_chunk,
                                             num_clusters_per_chunk))
        start = i * num_clusters_per_chunk
        end = start + num_clusters_in_chunk_i

        submitted = []
        # Create quiver bins and submit jobs
        all_todo = self.create_quiver_bins_and_submit_jobs(d=d, uc=uc,
                                                           partial_uc=partial_uc, refs=refs, keys=keys, start=start,
                                                           end=end, submitted=submitted, sge_opts=self.sge_opts)

        # Write submitted quiver jobs to
        # $root_dir/log/submitted_quiver_jobs.{i}of{num_chunks}.txt
        log_name = self.submitted_quiver_jobs_log_of_chunk_i(
            i=i, num_chunks=num_chunks)
        self.add_log("Writing submitted quiver jobs to {f}".format(f=log_name))
        with open(log_name, 'w') as f:
            f.write("\n".join(str(x[0]) + '\t' + str(x[1]) for x in submitted))

        # Write all quiver jobs of this workload to
        # $root_dir/log/quiver_jobs.{i}of{num_chunks}.sh
        sh_name = self.quiver_jobs_sh_of_chunk_i(i=i, num_chunks=num_chunks)
        self.add_log("Writing all quiver jobs to {f}".format(f=sh_name))
        with open(sh_name, 'w') as f:
            f.write("\n".join(["bash " + str(x) for x in all_todo]))

    def run(self):
        """Run quiver to polish all consensus isoforms predicted by ICE."""
        # Validate inputs
        self.validate_inputs()

        # One workload in total
        self.process_chunk_i(i=0, num_chunks=1)

        # Copy $root_dir/log/submitted_quiver_jobs.0of1.txt
        # to $root_dir/log/submitted_quiver_jobs.txt
        src = self.submitted_quiver_jobs_log_of_chunk_i(i=0, num_chunks=1)
        shutil.copyfile(src=src, dst=self.submitted_quiver_jobs_log)

        self.close_log()
        return 0


def add_ice_quiver_arguments(parser):
    """Add arguments for IceQuiver, not including IceQuiverPostprocess."""
    parser = add_cluster_root_dir_as_positional_argument(parser)
    parser = add_fofn_arguments(parser, bas_fofn=True, fasta_fofn=True)
    parser = add_sge_arguments(parser, quiver_nproc=True, blasr_nproc=True)
    return parser


# import os
# import sys
# from pbcore.util.ToolRunner import PBToolRunner
# from pbtools.pbtranscript.__init__ import get_version
#
#
# class IceQuiverRunner(PBToolRunner):
#
#     """IceQuiver runner."""
#
#     def __init__(self):
#         PBToolRunner.__init__(self, IceQuiver.desc)
#         add_ice_quiver_arguments(self.parser)
#
#     def getVersion(self):
#         """Get version string."""
#         return get_version()
#
#     def run(self):
#         """Run"""
#         logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
#                                                 v=self.getVersion()))
#         args = self.args
#         from pbtools.pbtranscript.ClusterOptions import SgeOptions
#         sge_opts = SgeOptions(os.getpid(),
#                               use_sge=args.use_sge,
#                               max_sge_jobs=args.max_sge_jobs,
#                               quiver_nproc=args.quiver_nproc)
#         try:
#             obj = IceQuiver(root_dir=args.root_dir,
#                             bas_fofn=args.bas_fofn,
#                             fasta_fofn=args.fasta_fofn,
#                             sge_opts=sge_opts)
#             obj.run()
#         except:
#             logging.exception("Exiting IceQuiver with return code 1.")
#             return 1
#         return 0
#
#
# def main():
#     """Main function"""
#     runner = IceQuiverRunner()
#     return runner.start()
#
# if __name__ == "__main__":
#     sys.exit(main())
