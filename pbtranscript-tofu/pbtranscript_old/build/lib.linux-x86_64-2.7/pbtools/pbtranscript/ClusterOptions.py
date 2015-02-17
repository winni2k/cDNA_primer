#!/usr/bin/env python
"""Define options for SGE configuration and the ICE algorithm."""


class SgeOptions(object):
    """Define options to configure SGE."""
    def __init__(self, unique_id, use_sge=False, max_sge_jobs=40,
                 blasr_nproc=24, gcon_nproc=8, quiver_nproc=8):
        self.unique_id = unique_id
        self.use_sge = use_sge
        self.max_sge_jobs = max_sge_jobs
        self.blasr_nproc = blasr_nproc
        self.gcon_nproc = gcon_nproc
        self.quiver_nproc = quiver_nproc

    def __str__(self):
        return "unqiueID={i}\n".format(i=self.unique_id) + \
               "use_sge={u}\n".format(u=self.use_sge) + \
               "max_sge_jobs={j}\n".format(j=self.max_sge_jobs) + \
               "blasr_nproc={n}\n".format(n=self.blasr_nproc) + \
               "gcon_nproc={n}\n".format(n=self.gcon_nproc) + \
               "quiver_nproc={t}\n".format(t=self.quiver_nproc)


class IceOptions(object):
    """Define ICE related options."""
    def __init__(self, cDNA_size="under1k", flnc_reads_per_split=20000,
        ece_penalty=1, ece_min_len=20, bestn=24, quiver=False,
        nfl_reads_per_split=30000):
        self.cDNA_size = str(cDNA_size)
        self.ece_penalty = int(ece_penalty)
        self.ece_min_len = int(ece_min_len)
        self.bestn = int(bestn)
        self.quiver = quiver
        # flnc reads per split
        self.flnc_reads_per_split = int(flnc_reads_per_split)
        # nfl reads per split
        self.nfl_reads_per_split = int(nfl_reads_per_split)

    @classmethod
    def cDNA_sizeBins(cls):
        """Return cDNA size bins."""
        return ("under1k", "between1k2k", "between2k3k", "above3k")

    @property
    def maxScore(self):
        """Return maximum blasr score according to estimated cDNA size."""
        if self.cDNA_size not in IceOptions.cDNA_sizeBins():
            raise ValueError("Invalid cDNA size: {cs}".
                format(cs=self.cDNA_size))
        d = {"under1k":-1000, "between1k2k":-2000, "between2k3k":-3000,
             "above3k":-5000}
        return d[self.cDNA_size]

    def __str__(self):
        return "cDNA_size={sz}\n".format(sz=self.cDNA_size) + \
               "maxScore={ms}\n".format(ms=self.maxScore) + \
               "ece_penalty={ep}\n".format(ep=self.ece_penalty) + \
               "ece_min_len={eml}\n".format(eml=self.ece_min_len) + \
               "bestn={bsn}\n".format(bsn=self.bestn) + \
               "quiver={quiver}\n".format(quiver=self.quiver) + \
               "flnc_reads_per_split={n}\n".format(n=self.flnc_reads_per_split) + \
               "nfl_reads_per_split={n}\n".format(n=self.nfl_reads_per_split)



