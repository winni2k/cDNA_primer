from pbtools.pbtranscript.ice.basQV_parallel import basQVcacher
from libc.math cimport log
from pbcore.io.FastaIO import FastaReader

class ProbFromQV:
    def __init__(self, input_fofn, fasta_filename=None, prob_threshold=.1, window_size=3):
        #self.qver = basQVcacherSparse(prob_threshold, window_size)
        self.qver = basQVcacher()
        self.input_fofn = input_fofn
        self.seqids = []
        self.prob_threshold = prob_threshold
        self.window_size = window_size

        with open(self.input_fofn) as f:
            for line in f:
                self.qver.add_bash5(line.strip())

        if fasta_filename is not None:
            self.add_seqs_from_fasta(fasta_filename)

    def get_smoothed(self, qID, qvname, position=None):
        return self.qver.get_smoothed(qID, qvname, position)

    def get(self, qID, qvname, position=None):
        return self.qver.get(qID, qvname, position)

    def add_seqs_from_fasta(self, fasta_filename, smooth=True):
        with open(fasta_filename) as f:
            newids = [r.name.split()[0] for r in FastaReader(f)]
        self.add_ids_from_fasta(newids, smooth)

    def add_ids_from_fasta(self, newids, smooth=True):
        self.qver.precache(newids)
        self.seqids += newids
        self.qver.presmooth(newids, self.window_size)
        #self.qver.remove_unsmoothed()

    def remove_ids(self, ids):
        for id in ids:
            self.seqids.remove(id)
            del self.qver.qv[id]

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        """
        """
        prob_sub = self.qver.get(qID, 'SubstitutionQV')
        prob_ins = self.qver.get(qID, 'InsertionQV')
        prob_del = self.qver.get(qID, 'DeletionQV')
        #prob_mat = 1 - prob_sub - prob_ins - prob_del  # not really right, but should be ok for now...
        # sanity check...sometimes this can have < 0 prob, so assign it a small prob like 0.001
        return calc_aln_log_prob(prob_sub, prob_ins, prob_del, len(prob_del), fakecigar, qStart, qEnd)


cdef double calc_aln_log_prob(list prob_sub, list prob_ins, list prob_del, int n, list fakecigar, int qStart, int qEnd):
    cdef int i, cur_q_pos
    cdef double score, tmp, one_three

#    prob_mat = []
#    for i in xrange(n):
#        tmp = 1 - prob_sub[i] - prob_ins[i] - prob_del[i] # not really right, but should be ok for now...
#        if tmp <= 0:
#            prob_mat.append(.001)
#        else:
#            prob_mat.append(tmp)

    one_three = log(1/3.)
    cur_q_pos = qStart
    score = 0.
    for x in fakecigar:
        if x == 'M':
            tmp = 1 - prob_sub[cur_q_pos] - prob_ins[cur_q_pos] - prob_del[cur_q_pos]
            if tmp <= 0: tmp = 0.001            
            score += log(tmp)
            cur_q_pos += 1
        elif x == 'S':
            score += log(prob_sub[cur_q_pos]) + one_three
            cur_q_pos += 1
        elif x == 'I':
            score += log(prob_ins[cur_q_pos]) + one_three
            cur_q_pos += 1
        else: # x == 'D', don't advance qpos
            score += log(prob_del[cur_q_pos])
    assert cur_q_pos == qEnd
#    del prob_mat
    return score


class fakeQVer:
    """
    Used by ProbFromModel to support the fake .get and .getsmoothed
    """
    def get(self, x, y, z=None): return 0.
    def get_smoothed(self, x, y, z=None): return 0.

class ProbFromModel:
    def __init__(self, r_mis, r_ins, r_del):
        self.qver = fakeQVer()
        self.r_mis = r_mis
        self.r_ins = r_ins
        self.r_del = r_del
        self.r_mat = 1 - r_mis - r_ins - r_del
        assert self.r_mat > 0

    def add_seqs_from_fasta(self, fasta_filename, smooth=True):
        pass # dummy, nothing to do

    def remove_ids(self, ids):
        pass # dummy, nothing to do

    def get_qver_smoothed(self, *args):
        return 0.

    def get(self, x, y, z=None):
        return self.qver.get(x, y, z)

    def get_smoothed(self, x, y, z=None):
        return self.qver.get_smoothed(x,y,z)

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        prob_mat = log(self.r_mat)
        prob_sub = log(self.r_mis)
        prob_ins = log(self.r_ins)
        prob_del = log(self.r_del)

        score = 0
        for x in fakecigar:
            if x == 'M':
                score += prob_mat
            elif x == 'S':
                score += prob_sub
            elif x == 'I':
                score += prob_ins
            else: # x == 'D', don't advance qpos
                score += prob_del
        return score
