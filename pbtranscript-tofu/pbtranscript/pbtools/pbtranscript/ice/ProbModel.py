"""Define Probability Models, including ProbFromQV and ProbFromModel."""
from pbcore.io.FastaIO import FastaReader
from pbcore.io.FastqIO import FastqReader
from pbtools.pbtranscript.io.BasQV import basQVcacher, fastqQVcacher
import numpy as np

class aloha:
    def __init__(self):
        pass

class ProbFromQV:
    """
    Probability model constructed from FOFN files using
    quality values.
    """
    def __init__(self, input_fofn, fasta_filename=None,
            prob_threshold=.03, window_size=5):

        self.qver = basQVcacher()
        self.input_fofn = input_fofn
        self.seqids = []
        self.prob_threshold = prob_threshold
        self.window_size = window_size
        self.full_prob = None

        with open(self.input_fofn) as f:
            for line in f:
                self.qver.add_bash5(line.strip())

        if fasta_filename is not None:
            self.add_seqs_from_fasta(fasta_filename)

        #self.qver.presmooth(self.seqids, self.window_size)

    def get_smoothed(self, qID, qvname, position=None):
        """
        Get smoothed QV of read=qID, type=qvname, position=position.
        """
        return self.qver.get_smoothed(qID, qvname, position)

    def get(self, qID, qvname, position=None):
        """
        Get QV of read=qID, type=qvname, position=position.
        """
        return self.qver.get(qID, qvname, position)

    def add_seqs_from_fasta(self, fasta_filename, smooth=True):
        """Add sequence ids from a fasta file."""
        with FastaReader(fasta_filename) as reader:
            newids = [r.name.split()[0] for r in reader]
        #with open(fasta_filename) as f:
        #    newids = [r.id for r in SeqIO.parse(f, 'fasta')]
        self.add_ids_from_fasta(newids, smooth)

    def add_ids_from_fasta(self, newids, _smooth=True):
        """Add sequence ids."""
        self.qver.precache(newids)
        self.seqids += newids
        self.qver.presmooth(newids, self.window_size)
        #self.qver.remove_unsmoothed()

    def remove_ids(self, ids):
        """Remove ids from self.seqids."""
        for _id in ids:
            self.seqids.remove(_id)
            del self.qver.qv[_id]

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        """
        aln_record is a pysam.AlignedRead
        """
        # using sparse matrix must make everything positive
        prob_sub = self.qver.get(qID, 'SubstitutionQV') + .0001
        prob_ins = self.qver.get(qID, 'InsertionQV')  + .0001
        prob_del = self.qver.get(qID, 'DeletionQV') + .0001

        # not really right, but should be ok for now...
        prob_mat = 1 - prob_sub - prob_ins - prob_del
        # sanity check...sometimes this can have < 0 prob,
        # so assign it a small prob like 0.001
        ii = (prob_mat <= 0)
        prob_mat[ii.nonzero()] = .001
        prob_sub = np.log(prob_sub)
        prob_ins = np.log(prob_ins)
        prob_del = np.log(prob_del)
        prob_mat = np.log(prob_mat)

        one_three = np.log(1/3.)
        cur_q_pos = qStart
        score = 0
        for x in fakecigar:
            if x == 'M':
                score += prob_mat[cur_q_pos]
                cur_q_pos += 1
            elif x == 'S':
                score += prob_sub[cur_q_pos] + one_three
                cur_q_pos += 1
            elif x == 'I':
                score += prob_ins[cur_q_pos] + one_three
                cur_q_pos += 1
            else: # x == 'D', don't advance qpos
                score += prob_del[cur_q_pos]
        assert cur_q_pos == qEnd
        return score

class ProbFromFastq:
    """
    Probability model constructed from Fastq files using a single QV for everything
    """
    def __init__(self, fastq_filename,
            prob_threshold=.03, window_size=5):

        self.qver = fastqQVcacher()
        self.fastq_filename = fastq_filename
        self.seqids = []
        self.prob_threshold = prob_threshold
        self.window_size = window_size
        self.full_prob = None

        self.add_seqs_from_fastq(fastq_filename)

    def get_smoothed(self, qID, qvname, position=None):
        """
        Get smoothed QV of read=qID, type=qvname, position=position.
        In reality, qvname is ignored.
        """
        return self.qver.get_smoothed(qID, qvname, position)

    def get(self, qID, qvname, position=None):
        """
        Get QV of read=qID, type=qvname, position=position.
        In reality, qvname is ignored.
        """
        return self.qver.get(qID, qvname, position)

    def add_seqs_from_fastq(self, fastq_filename, smooth=True):
        """Add sequence ids from a fastq file."""
        self.qver.precache_fastq(fastq_filename)
        newids = [r.name.split()[0] for r in FastqReader(fastq_filename)]
        self.qver.presmooth(newids, self.window_size)
        #self.qver.remove_unsmoothed()

    def remove_ids(self, ids):
        """Remove ids from self.seqids."""
        for _id in ids:
            self.seqids.remove(_id)
            del self.qver.qv[_id]

    def calc_prob_from_aln(self, qID, qStart, qEnd, fakecigar):
        """
        aln_record is a pysam.AlignedRead
        """
        prob_err = self.qver.get(qID, None) + .0001

        # not really right, but should be ok for now...
        prob_mat = 1 - prob_err 
        # sanity check...sometimes this can have < 0 prob,
        # so assign it a small prob like 0.001
        ii = (prob_mat <= 0)
        prob_mat[ii.nonzero()] = .001
        prob_err = np.log(prob_err)
        prob_mat = np.log(prob_mat)

        one_three = np.log(1/3.)
        cur_q_pos = qStart
        score = 0
        for x in fakecigar:
            if x == 'M':
                score += prob_mat[cur_q_pos]
                cur_q_pos += 1
            elif x == 'S':
                score += prob_err[cur_q_pos] + one_three
                cur_q_pos += 1
            elif x == 'I':
                score += prob_err[cur_q_pos] + one_three
                cur_q_pos += 1
            else: # x == 'D', don't advance qpos
                score += prob_err[cur_q_pos]
        assert cur_q_pos == qEnd
        return score


class fakeQVer:
    """
    Used by ProbFromModel to support the fake .get and .getsmoothed
    """
    def __init__(self):
        pass

    def get(self, _x, _y, _z=None):
        """Fake get()."""
        return 0.

    def get_smoothed(self, qID, qvname, position=None):
        """Fake get_smoothed."""
        return 0.


class ProbFromModel:
    """Probability model from fixed indel/substitution rates."""
    def __init__(self, r_mis, r_ins, r_del):
        self.qver = fakeQVer()
        self.r_mis = r_mis
        self.r_ins = r_ins
        self.r_del = r_del
        self.r_mat = 1 - r_mis - r_ins - r_del
        assert self.r_mat > 0

    def add_seqs_from_fasta(self, fasta_filename, smooth=True):
        """# dummy, nothing to do"""
        pass

    def remove_ids(self, ids):
        """# dummy, nothing to do"""
        pass

    #def get_qver_smoothed(self, *args):
    #    return 0.

    def get(self, x, y, z=None):
        """Return QV for x, y, z."""
        return self.qver.get(x, y, z)

    def get_smoothed(self, qID, qvname, position=None):
        """Get smoothed QV for qId, qvname, position."""
        return self.qver.get_smoothed(qID=qID, qvname=qvname,
                                      position=position)

    def calc_prob_from_aln(self, _qID, _qStart, _qEnd, fakecigar):
        """Calculate probability from an alingment."""
        prob_mat = np.log(self.r_mat)
        prob_sub = np.log(self.r_mis)
        prob_ins = np.log(self.r_ins)
        prob_del = np.log(self.r_del)

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



