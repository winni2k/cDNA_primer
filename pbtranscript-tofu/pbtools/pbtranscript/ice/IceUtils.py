"""Define functions useful for IceInit and IceIterative."""
import os
import os.path as op
import glob
import logging
import shutil
import filecmp
import random
import numpy as np
from pbcore.util.Process import backticks
from pbcore.io import FastaReader, FastaWriter
from pbtools.pbtranscript.Utils import realpath, mkdir, \
    get_files_from_fofn, write_files_to_fofn
from pbtools.pbtranscript.io.BLASRRecord import BLASRM5Reader
from pbtools.pbtranscript.findECE import findECE

__author__ = 'etseng@pacificbiosciences.com'

#define gcon script for ice.
gcon_py = "pbdagcon_wrapper.py"

# Define data sets for sge sanity check.
dataDir = op.join(op.dirname(op.dirname(op.realpath(__file__))), "data")
GCON_IN_FA = op.join(dataDir, "gcon_in.fa")
GCON_OUT_FA = op.join(dataDir, "gcon_out.fa")

def clean_up_after_ICE(dirname):
    """
    Clean up files that are not needed after ICE is successfully completed    
    """
    patterns = ['*.blasr', 'ref_consensus*', 'current.fasta*', 'input.split*fa']
    for pattern in patterns:
        for file in glob.iglob(os.path.join(dirname, pattern)):
            try:
                os.remove(file)
            except:
                logging.error("Unable to clean up file {0}.".format(file))
        

def sanity_check_gcon():
    """Sanity check gcon."""
    cmd = gcon_py + " --help"
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        msg = gcon_py + " is not installed."
        raise RuntimeError(msg)
    return gcon_py

def sanity_check_sge(scriptDir, testDirName="gcon_test_dir"):
    """Sanity check if sge can work."""
    scriptDir = realpath(scriptDir)
    testDir = op.join(scriptDir, testDirName)

    if not op.exists(scriptDir):
        os.makedirs(scriptDir)
    if not op.exists(testDir):
        os.makedirs(testDir)

    testSh = op.join(scriptDir, 'test.sh')
    consensusFa = op.join(testDir, "g_consensus.fa")
    testInFa = op.join(testDir, "gcon_in.fa")
    shutil.copy(GCON_IN_FA, testInFa)
    assert(op.exists(testInFa))

    with open(testSh, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("{gcon}".format(gcon=gcon_py) +
                " {inFa}".format(inFa=testInFa) +
                " {testDir}/g_consensus".format(testDir=testDir) +
                " c1\n")

    assert(op.exists(testSh))
    cmd = "qsub -sync y -pe smp 1 -cwd -S /bin/bash -V " + \
          "-e /dev/null -o /dev/null {t}".format(t=testSh)
    logging.debug("Submitting cmd: " + cmd)
    _out, _code, _msg = backticks(cmd)
    
    answer = FastaReader(GCON_OUT_FA).__iter__().next()
    tester = FastaReader(consensusFa).__iter__().next()

    if answer.name!=tester.name or answer.sequence!=tester.sequence:
        errMsg = "Trouble running qsub or output is not as " + \
                 "expected ({0} and {1} must agree). Abort!".format(
                 consensusFa, GCON_OUT_FA)
        logging.error(errMsg)
        return False
    else:
        shutil.rmtree(testDir)
        logging.info("sge and gcon check passed.")
        return True


def eval_blasr_alignment(record, qver_get_func,
        sID_starts_with_c, qv_prob_threshold):
    """
    Takes a BLASRRecord (blasr -m 5) and goes through the
    alignment string
    ex: |||**||||**|||*|*|
    to determine the sequence of 'M' (matches), 'S' (sub), 'I', 'D'

    qver_get_func --- could be either basQV.basQVcacher.get() or
        basQV.basQVcacher.get_smoothed()

    For any non-match, if either or both query/target's QV indicate
    that the event ('S', 'I', 'D') is expected
    (ex: insertion prob >= qv_prob_threshold),
    then it does not count as a penalty.

    Returns: cigar string, binary ECE array

    NOTE: long insertions/deletions are still a difficult problem
    because alignments can be arbitrary right now the quick solution is:
          use probqv get_smoothed
    however, with homopolymers, penalization can still happen unless
    I write code to check specifically for homopolymers, (otherwise the
    cigar_str[-1]=='D' or 'I' sets in). -- Liz
    """
    q_index = 0
    s_index = 0
    cigar_str = ''
    # binary array of 0|1 where 1 is a penalty
    ece = np.zeros(len(record.alnStr), dtype=np.int)
    #pdb.set_trace()
    for offset, nt_aln in enumerate(record.alnStr):
        if nt_aln == '|': # match
            cigar_str += 'M'
            q_index += 1
            s_index += 1
        elif record.qAln[offset] == '-': # deletion
            if cigar_str[-1] == 'D' or (\
                qver_get_func(record.qID, 'DeletionQV', q_index+1) < qv_prob_threshold and \
               (sID_starts_with_c or qver_get_func(record.sID, 'InsertionQV', s_index) < qv_prob_threshold)):
                # case 1: last one was also a D (so q_index did not advance)
                # case 2: both QVs were good yet still a non-match, penalty!
                ece[offset] = 1
            cigar_str += 'D'
            s_index += 1
        elif record.sAln[offset] == '-': # insertion
            if cigar_str[-1] == 'I' or (\
                qver_get_func(record.qID, 'InsertionQV', q_index) < qv_prob_threshold and \
                (sID_starts_with_c or qver_get_func(record.sID, 'DeletionQV', s_index+1) < qv_prob_threshold)):
                # case 1: last one was also a I (so s_index did not advance)
                # case 2: both QVs were good yet still a no-match
                ece[offset] = 1
            cigar_str += 'I'
            q_index += 1
        else: # substitution
            cigar_str += 'S'
            if qver_get_func(record.qID, 'SubstitutionQV', q_index) < qv_prob_threshold and \
               (sID_starts_with_c or qver_get_func(record.sID, 'SubstitutionQV', s_index) < qv_prob_threshold):
                ece[offset] = 1
            q_index += 1
            s_index += 1
    #assert q_index == len(record.qAln) - record.qAln.count('-')
    #assert s_index == len(record.sAln) - record.sAln.count('-')
    return cigar_str, ece


class HitItem(object):
    """
    Simply define an object class for saving items produced by
    blasr_against_ref.
    """
    def __init__(self, qID, cID, qStart=None, qEnd=None,
            missed_q=None, missed_t=None,
            fakecigar=None, ece_arr=None):
        self.qID = qID
        self.cID = cID
        self.qStart = qStart
        self.qEnd = qEnd
        self.missed_q = missed_q
        self.missed_t = missed_t
        self.fakecigar = fakecigar
        self.ece_arr = ece_arr


def blasr_against_ref(output_filename, is_FL, sID_starts_with_c,
    qver_get_func, qv_prob_threshold=.1,
    ece_penalty=1, ece_min_len=20, same_strand_only=True):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    (3) less than 90% aligned or more than 50 bp missed

    qver_get_func --- should be basQV.basQVcacher.get() or
                      .get_smoothed(), or can just pass in
                      lambda (x, y): 1. to ignore QV
    """
    with BLASRM5Reader(output_filename) as reader:
        for r in reader:
            missed_q = r.qStart + r.qLength - r.qEnd
            missed_t = r.sStart + r.sLength - r.sEnd

            if sID_starts_with_c:
                # because all consensus should start with
                # c<cluster_index>
                assert r.sID.startswith('c')
                if r.sID.find('/') > 0:
                    r.sID = r.sID.split('/')[0]
                if r.sID.endswith('_ref'):
                    # probably c<cid>_ref
                    cID = int(r.sID[1:-4])
                else:
                    cID = int(r.sID[1:])
            else:
                cID = r.sID

            # self hit, useless!
            # low identity not allowed
            # opposite strand not allowed!
            if (cID == r.qID or
                r.identity < 90. or
                (r.strand == '-' and same_strand_only)):
                yield HitItem(qID=r.qID, cID=cID)
                continue

            # full-length case: allow up to 100bp of 5' not aligned
            # and 30bp of 3' not aligned
            # non-full-length case: not really tested...don't use
            if is_FL and (r.sStart > 100 or r.qStart > 100 or
                          (r.sLength-r.sEnd > 30) or
                          (r.qLength-r.qEnd > 30)):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                cigar_str, ece_arr = eval_blasr_alignment(
                    record=r,
                    qver_get_func=qver_get_func,
                    sID_starts_with_c=sID_starts_with_c,
                    qv_prob_threshold=qv_prob_threshold)

                if alignment_has_large_nonmatch(ece_arr,
                    ece_penalty, ece_min_len):
                    yield HitItem(qID=r.qID, cID=cID)
                else:
                    yield HitItem(qID=r.qID, cID=cID,
                        qStart=r.qStart, qEnd=r.qEnd,
                        missed_q=missed_q*1./r.qLength,
                        missed_t=missed_t*1./r.sLength,
                        fakecigar=cigar_str,
                        ece_arr=ece_arr)


def alignment_has_large_nonmatch(ece_arr, penalty, min_len):
    """
    penalty of (-)1: 50%
    penalty of (-)2: 66%
    penalty of (-)4: 80%
    penalty of (-)9: 90%

    Return True when alignment has large non-matches not explained
    by low base QVs (in other words, "reject" as an isoform hit and
    don't put in the same cluster)
    """
    ece_arr = ece_arr * (penalty + 1)
    s = [0] + list(ece_arr - penalty)
    # fix this later to something faster & better
    return (len(findECE(s, len(s), min_len, True)) > 0)


def possible_merge(r, ece_penalty, ece_min_len):
    """
    r --- BLASRM5Record
    Criteria:
    (1) identity >= 99% and same strand
    (2) check criteria for how much is allowed to differ on the
        5' / 3' ends

    NOTE: 100 bp on 5' and 30 bp on 3' leniency is now HARD-CODED!
    Should change later
    """
    if r.sID == r.qID or r.identity < 99 or r.strand == '-':
        return False
    # intentional here to prevent disrupting future ICE runs
    # MORE lenient on 5' but NOT on 3'
    if ((r.qLength - r.qEnd) > 30 or (r.sLength - r.sEnd) > 30 or \
        r.qStart > 100 or r.sStart > 100):
        return False

    arr = np.array([(x == '*') * 1 for x in r.alnStr])
    if alignment_has_large_nonmatch(ece_arr=arr,
                                    penalty=ece_penalty,
                                    min_len=ece_min_len):
        return False
    return True


def get_the_only_fasta_record(fa):
    """Input fasta file should contain exactly one FastaRecord,
    return the fastas record."""
    rs = [r for r in FastaReader(fa)]
    if len(rs) != 1:
        errMsg = "Cluster fasta file {fa} must contain only one read.".\
            format(fa=fa)
        raise ValueError(errMsg)
    return rs[0]


"""
The following methods was originally created by jchin:
    /home/UNIXHOME/jchin/depot_mp27/jchin/rset_quvier.py
, and then modified by etseng.

Input: input.fasta.fofn (shared),
       per-cluster in.fa,
       per-cluster g_consensus.fa

-- input.fasta.fofn should be raw fasta files
   (pls2fasta -maskRegion) of input.fofn

Within each cluster:
1) create in.raw.fa based on input.fasta.fofn & in.fa,
   putting in raw (unrolled) fasta of each ZMW
2) blasr (1) to g_consensus.fa output as SAM

This is faster than using regions.fofn because it still reads
through the whole .bax.h5 files
"""
def is_blank_sam(samfile):
    """
    return True if the SAM file only has @xx header and NO alignment
    """
    with open(samfile) as f:
        for line in f:
            if not line.startswith('@'):
                return False
    return True

def concat_sam(samfiles, outsam_filename):
    """
    Header looks like:
    @HD     VN:1.3.1
    @SQ     SN:c31  LN:3104 M5:ef7d3f84dea9d9face43e6fd5b6336c4
    @RG     ID:2caa54eef6   PU:in.raw_with_partial.fa       SM:NO_CHIP_ID
    @PG     ID:BLASR        VN:1.3.1.126469 CL:blasr in.raw_with_partial.fa g_consensus.fa -nproc 12 -bestn 5 -nCandidates 10 -sam -out out.sam

    NOTE: check for M5 conflicts; manipulate them if it conflicts
    """
    f_sq = open(outsam_filename + '.sq', 'w')
    f_bd = open(outsam_filename + '.bd', 'w')

    rg_line = None
    pg_line = None

    md5_seen = set()

    h = open(samfiles[0])
    line = h.readline()
    assert line.startswith('@HD')
    f_sq.write(line)
    line = h.readline()
    assert line.startswith('@SQ')
    line = h.readline()
    assert line.startswith('@RG')
    rg_line = line # write at the end
    line = h.readline()
    assert line.startswith('@PG')
    pg_line = line # write at the end
    h.close()

    for f in samfiles:
        with open(f) as h:
            assert h.readline().startswith('@HD')
            line = h.readline()
            assert line.startswith('@SQ')
            # ------- check for MD5 conflicts ----------- #
            m5 = line.strip().split()[-1]
            assert m5.startswith("M5:")
            if m5 not in md5_seen:
                f_sq.write(line)
                md5_seen.add(m5)
            else:
                s = list(m5[3:])
                while True:
                    # create a random m5 string.
                    random.shuffle(s)
                    s = "".join(s)
                    if s not in md5_seen:
                        break
                line = line[:line.find('M5:')] + 'M5:' + s + '\n'
                logging.debug("MD5 conflict: change to {0}".format(s))
                md5_seen.add(s)
                f_sq.write(line)
            # ----- end MD5 checking and writing --------- #
            assert h.readline().startswith('@RG')
            assert h.readline().startswith('@PG')
            for line in h:
                f_bd.write(line)

    f_bd.close()
    f_sq.write(rg_line)
    f_sq.write(pg_line)
    f_sq.close()

    cmd = "cat {0}.sq {0}.bd > {0}".format(outsam_filename)
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        raise IOError("Failed to concat sam files! Abort." + _msg)

    os.remove(f_sq.name)
    os.remove(f_bd.name)


def convert_fofn_to_fasta(fofn_filename, out_filename, fasta_out_dir,
                          force_overwrite=False):
    """
    For each .bax.h5 file, create .bax.h5.fasta file and save paths to
    out_filename, which should usually be 'input.fasta.fofn'
    """
    logging.info("Converting fofn {fofn} to fasta.".format(fofn=fofn_filename))
    in_fns = get_files_from_fofn(fofn_filename)
    out_fns = []
    mkdir(fasta_out_dir)
    for in_fn in in_fns:
        logging.debug("converting h5 file: {f}.".format(f=in_fn))
        if not (in_fn.endswith('.bax.h5') or in_fn.endswith('.bas.h5')):
            raise ValueError("fofn file {fofn} ".format(fofn=fofn_filename) +
                             "should only contain bax/bas.h5 files.")

        # e.g. m111xxxx.1.bax.h5 ==>
        #      tmp_out_file = m11xxxx.1.bax.h5.fasta.tmp
        #      out_file = m11xxxx.1.bax.h5.fasta
        in_basename = op.basename(in_fn)
        tmp_out_file = op.join(fasta_out_dir, in_basename + '.fasta.tmp')
        out_file = op.join(fasta_out_dir, in_basename + '.fasta')
        if op.exists(out_file) and not force_overwrite:
            logging.debug("File {0} already exists. skipping.".format(out_file))
        else:
            cmd = "pls2fasta {0} {1} ".format(in_fn, tmp_out_file) + \
                  "-minSubreadLength 300 -minReadScore 750 -trimByRegion"
            logging.debug("CMD: {cmd}".format(cmd=cmd))
            _out, _code, _msg = backticks(cmd)
            if _code != 0:
                raise RuntimeError("CMD failed: {cmd}\n".format(cmd=cmd) + _msg)
            trim_subread_flanks(tmp_out_file, out_file)
        out_fns.append(out_file)
        if op.exists(tmp_out_file):
            os.remove(tmp_out_file)
    write_files_to_fofn(out_fns, out_filename)


def trim_subread_flanks(fasta_filename, output_filename,
                        trim_len=100, min_len=100):
    """
    fasta_filename --- should be subread output from pls2fasta

    trim first/last 100bp (which contains primer&polyA) away and correct
    coordinates
    """
    with FastaWriter(output_filename) as writer, \
         FastaReader(fasta_filename) as reader:
        for r in reader:
            # ex: m14011..._s1_p0/15/1305_4354
            movie, hn, s_e = r.name.split()[0].split('/')
            s, e = s_e.split('_')
            s, e = int(s), int(e)
            assert s < e
            s2 = s + trim_len
            e2 = e - trim_len
            if e2 - s2 >= min_len:
                newname = "{0}/{1}/{2}_{3}".format(movie, hn, s2, e2)
                newseq = r.sequence[trim_len:-trim_len]
                writer.writeRecord(newname, newseq)

def build_sa(input_fasta, out_sa):
    """Generate suffix array of input_fasta"""
    if op.exists(input_fasta):
        cmd = "sawriter {o} {i} -blt 8 -welter ".format(o=out_sa, i=input_fasta)
        _out, _code, _msg = backticks(cmd)
        if _code == 0:
            return True
        else:
            # If failed to generate suffix array, warning.
            logging.warn("Unable to create suffix array for {f}.".format(f=input_fasta))
            return False
    else:
        raise IOError("Unable to find fasta file {f}.".format(f=input_fasta))


def write_in_raw_fasta(input_fasta_d, in_seqids,
                       out_fa, ignore_keyerror=False):
    """
    input_fasta_d --- miscBio.MetaFastaReader
    input fasta should be in format <movie>/<holeNumber>/<subread or CCS stuff>

    Create a out_fa fasta where we dump the "raw" (unrolled) of every ZMW from in_fa
    """
    movies = set()
    zmw_seen = set()
    with open(out_fa, 'w') as f:
        for seqid in in_seqids:
            try:
                zmw = seqid[:seqid.rfind('/')]
                if zmw not in zmw_seen:
                    movies.add(zmw.split('/')[0])
                    for rec in input_fasta_d[zmw]:
                        f.write(">{0}\n{1}\n".format(rec.name, rec.sequence))
                    #f.write(">{0}\n{1}\n".format(zmw, input_fasta_d[zmw].seq))
                    zmw_seen.add(zmw)
            except KeyError:
                if ignore_keyerror:
                    logging.debug("Ignorning non-existent key " + zmw)
                else:
                    raise ValueError, "{0} doesn't exist. Abort!".format(zmw)
    return movies


def blasr_sam_for_quiver(input_fasta, ref_fasta,
                         out_sam_filename,
                         run_cmd=True, blasr_nproc=12):
    """
    #input_fofn --- should be input.fofn
    input_fasta --- should be in.raw.fa
    ref_fasta --- reference fasta (ex: g_consensus.fa) to align to
    #output_dir --- if None, automatically set to where ref_fasta is

    run blasr -clipping soft to get sam
    """
    #if output_dir is None:
    #    output_dir = op.dirname(ref_fasta)
    #if movies is not None:
    #    f = open(input_fasta + '.fofn', 'w')
    #    for line in open(input_fofn):
    #        if op.basename(line).split('.')[0] in movies:
    #            f.write(line)
    #    f.close()
    #    input_fofn = f.name
    #out_sam = op.join(output_dir, out_sam_filename)
    #TODO: review code

    cmd = "blasr {i} ".format(i=input_fasta) + \
          "{r} ".format(r=ref_fasta) + \
          "-nproc {n} ".format(n=blasr_nproc) + \
          "-bestn 5 -nCandidates 10 -sam -clipping soft " + \
          "-out {o}".format(o=out_sam_filename)
    logging.debug("CMD: " + cmd)
    if run_cmd:
        _out, _code, _msg = backticks(cmd)
        if _code != 0:
            raise RuntimeError("CMD failed: {cmd}\n{e}".
                format(cmd=cmd, e=_msg))
    return cmd


