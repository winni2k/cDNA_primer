#!/usr/bin/env python
"""Define options for SGE configuration and the ICE algorithm."""


class SgeOptions(object):

    """Define options to configure SGE."""

    def __init__(self, unique_id, use_sge=False, max_sge_jobs=40,
                 blasr_nproc=24, gcon_nproc=8, quiver_nproc=8, sge_env_name='smp', sge_queue=None):
        self.unique_id = unique_id
        self.use_sge = use_sge
        self.max_sge_jobs = max_sge_jobs
        self.blasr_nproc = blasr_nproc
        self.gcon_nproc = gcon_nproc
        self.quiver_nproc = quiver_nproc
        self.sge_env_name = sge_env_name
        self.sge_queue = sge_queue

    def __str__(self):
        return "unqiueID={i}\n".format(i=self.unique_id) + \
               "use_sge={u}\n".format(u=self.use_sge) + \
               "max_sge_jobs={j}\n".format(j=self.max_sge_jobs) + \
               "blasr_nproc={n}\n".format(n=self.blasr_nproc) + \
               "gcon_nproc={n}\n".format(n=self.gcon_nproc) + \
               "quiver_nproc={t}\n".format(t=self.quiver_nproc) + \
               "sge_env_name={0}\n".format(self.sge_env_name) + \
               "sge_queue={0}\n".format(self.sge_queue)

    def cmd_str(self, show_blasr_nproc=False, show_gcon_nproc=False,
                show_quiver_nproc=False, show_sge_env_name=False, show_sge_queue=False):
        """Return a cmd string."""
        cmd = ""
        if self.use_sge is True:
            cmd += "--use_sge "
            cmd += "--max_sge_jobs={n} ".format(n=self.max_sge_jobs)
            cmd += "--unique_id={n} ".format(n=self.unique_id)
        if show_blasr_nproc is True and self.blasr_nproc is not None:
            cmd += "--blasr_nproc={n} ".format(n=self.blasr_nproc)
        if show_gcon_nproc is True and self.gcon_nproc is not None:
            cmd += "--gcon_nproc={n} ".format(n=self.gcon_nproc)
        if show_quiver_nproc is True and self.quiver_nproc is not None:
            cmd += "--quiver_nproc={n} ".format(n=self.quiver_nproc)
        if show_sge_env_name:
            cmd += "--sge_env_name={0}".format(self.sge_env_name)
        if show_sge_queue and self.sge_queue is not None:
            cmd += "--sge_queue={0}".format(self.sge_queue)
        return cmd


class IceOptions(object):

    """Define ICE related options."""

    def __init__(self, flnc_reads_per_split=20000,
            ece_penalty=1, ece_min_len=20, bestn=24, quiver=False, use_finer_qv=False, targeted_isoseq=False):
                 #nfl_reads_per_split=30000):
        self.cDNA_size = "under1k"   # set a default, this will be automatically changed in IceIterative later
        self.ece_penalty = int(ece_penalty)
        self.ece_min_len = int(ece_min_len)
        self.bestn = int(bestn)
        self.quiver = quiver
        # flnc reads per split
        self.flnc_reads_per_split = int(flnc_reads_per_split)
        self.use_finer_qv = use_finer_qv
        # nfl reads per split (Liz: moving this to group with --nfl_fa)
        #self.nfl_reads_per_split = int(nfl_reads_per_split)
        self.targeted_isoseq = targeted_isoseq

    @classmethod
    def cDNA_sizeBins(cls):
        """Return cDNA size bins."""
        return ("under1k", "between1k2k", "between2k3k", "3to5k", "above5k")

    @property
    def maxScore(self):
        """Return maximum blasr score according to estimated cDNA size."""
        if self.cDNA_size not in IceOptions.cDNA_sizeBins():
            raise ValueError("Invalid cDNA size: {cs}".
                             format(cs=self.cDNA_size))
        d = {"under1k": -1000, "between1k2k": -2000, "between2k3k": -3000,
             "3to5k": -5000, "above5k": -8000}
        return d[self.cDNA_size]

    @property
    def minLength(self):
        """
        Used by DALIGNER for minimum length
        """
        if self.cDNA_size not in IceOptions.cDNA_sizeBins():
            raise ValueError("Invalid cDNA size: {cs}".
                             format(cs=self.cDNA_size))
        d = {"under1k": 300, "between1k2k": 800, "between2k3k": 1500,
             "3to5k": 2500, "above5k": 3500}
        return d[self.cDNA_size]

    def __str__(self):
        return "cDNA_size={sz}\n".format(sz=self.cDNA_size) + \
               "maxScore={ms}\n".format(ms=self.maxScore) + \
               "ece_penalty={ep}\n".format(ep=self.ece_penalty) + \
               "ece_min_len={eml}\n".format(eml=self.ece_min_len) + \
               "bestn={bsn}\n".format(bsn=self.bestn) + \
               "quiver={quiver}\n".format(quiver=self.quiver) + \
               "flnc_reads_per_split={n}\n".\
            format(n=self.flnc_reads_per_split) + \
               "use_finer_qv={qv}\n".format(qv=self.use_finer_qv)# + \
               #"nfl_reads_per_split={n}\n".format(n=self.nfl_reads_per_split)


class IceQuiverHQLQOptions(object):

    """Define HQ/LQ isoforms related options"""

    def __init__(self, qv_trim_5=200, qv_trim_3=50, hq_quiver_min_accuracy=0.99,
                 hq_isoforms_fa=None, hq_isoforms_fq=None,
                 lq_isoforms_fa=None, lq_isoforms_fq=None):
        # Ignore QV of n bases in the 5' end
        self.qv_trim_5 = int(qv_trim_5)
        # Ignore QV of n bases in the 3' end
        self.qv_trim_3 = int(qv_trim_3)
        # Minimum allowed quiver accuracy to mark an isoform as HQ
        self.hq_quiver_min_accuracy = float(hq_quiver_min_accuracy)

        self.hq_isoforms_fa = hq_isoforms_fa
        self.hq_isoforms_fq = hq_isoforms_fq
        self.lq_isoforms_fa = lq_isoforms_fa
        self.lq_isoforms_fq = lq_isoforms_fq

    def __str__(self):
        return "qv_trim_5={n}\n".format(n=self.qv_trim_5) + \
               "qv_trim_3={n}\n".format(n=self.qv_trim_3) + \
               "hq_quiver_min_accuracy={n}\n".\
            format(n=self.hq_quiver_min_accuracy) + \
               "HQ isoforms fasta={fa}\n".format(fa=self.hq_isoforms_fa) + \
               "HQ isoforms fastq={fq}\n".format(fq=self.hq_isoforms_fq) + \
               "LQ isoforms fasta={fa}\n".format(fa=self.lq_isoforms_fa) + \
               "LQ isoforms fastq={fq}\n".format(fq=self.lq_isoforms_fq)

    def cmd_str(self):
        """Return a cmd string."""
        cmd = "--hq_quiver_min_accuracy={n} ".\
                format(n=self.hq_quiver_min_accuracy) + \
              "--qv_trim_5={n} ".format(n=self.qv_trim_5) + \
              "--qv_trim_3={n} ".format(n=self.qv_trim_3)
        if self.hq_isoforms_fa is not None:
            cmd += "--hq_isoforms_fa={f} ".format(f=self.hq_isoforms_fa)
        if self.hq_isoforms_fq is not None:
            cmd += "--hq_isoforms_fq={f} ".format(f=self.hq_isoforms_fq)
        if self.lq_isoforms_fa is not None:
            cmd += "--lq_isoforms_fa={f} ".format(f=self.lq_isoforms_fa)
        if self.lq_isoforms_fq is not None:
            cmd += "--lq_isoforms_fq={f} ".format(f=self.lq_isoforms_fq)
        return cmd

