import os, sys, random, subprocess
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from pbcore.io import BasH5Reader
import GFF

def select_random_sequences(bas_filename, out_fa, out_fq, random_prob, num_of_seqs, use_CCS=False, min_seq_len=500):
    """
    random_prob --- prob between (0, 1] to use for selecting sequences
    num_of_seqs --- selection stops either by end of .bas.h5 or reaches this limit
    use_CCS --- if True, selects CCS; other selects subreads
    """
    bas = BasH5Reader(bas_filename)
    count = 0
    for zmw_number in bas.sequencingZmws:
        zmw = bas[zmw_number]
        if use_CCS:
            if zmw.ccsRead is not None and len(zmw.ccsRead.basecalls()) >= min_seq_len and random.random() <= random_prob:
                s = zmw.ccsRead
                count += 1
                out_fa.write(">{0}\n{1}\n".format(s.readName, s.basecalls()))
                out_fq.write("@{0}\n{1}\n+\n{2}\n".format(s.readName, s.basecalls(),  "".join(chr(x+33) for x in s.QualityValue())))
        else:
            for subread in zmw.subreads:
                if len(subread.basecalls()) >= min_seq_len and random.random() <= random_prob:
                    s = subread
                    count += 1
                    out_fa.write(">{0}\n{1}\n".format(s.readName, s.basecalls()))
                    out_fq.write("@{0}\n{1}\n+\n{2}\n".format(s.readName, s.basecalls(),  "".join(chr(x+33) for x in s.QualityValue())))
        if count >= num_of_seqs:
            return count
    return count



def calc_overlap(r1, r2):
    """
    Used for seeing whether some of the GMAP chimeras are actually missed adapters
    Criteria:
    (1) same chromosome
    (2) overlap is at least 50% of both records
    """
    if r1.chr != r2.chr: return False
    a = max(r1.start, r2.start)
    b = min(r1.end, r2.end)
    if a < b and (b-a)>=(r1.end-r1.start)*.5 and (b-a)>=(r2.end-r2.start)*.5: return True
    return False


def run_gmap(fa_filename, fq_filename=None, gmap_filename=None, tempdir='/scratch'):
    """
    Run gmap on FASTA file, extracting the alignment identity and comparing it with FASTQ QVs
    """
    if os.path.exists(gmap_filename):
        print >> sys.stderr, gmap_filename, "exists. No need to run gmap again."
        out_filename = gmap_filename
    else:
        out_filename = os.tempnam(tempdir)
        cmd = "gmap -D /home/UNIXHOME/etseng/share/gmap_db -d {db} -t 12 -f gff3_gene -n 0 {i} > {o}\n".format(i=fa_filename, o=out_filename, db=DBNAME)
        if subprocess.check_call(cmd, shell=True) !=0 :
            print >> sys.stderr, "Error running: {0}. Abort.".format(cmd)
            sys.exit(-1)
        
    tally = defaultdict(lambda: []) # perZMW --> gmapRecord
    coverages = defaultdict(lambda: 0)
    for r in GFF.gmapGFFReader(out_filename):
        tally[r.seqid].append(r)
    
    obs2exp = [] # (exp. avg. accuracy, obs. avg. accuracy, size)
    unmapped = 0 # count of not aligned by GMAP
    coverages = defaultdict(lambda: 0)
    zmw_seen = set()
    if fq_filename is None:
        for r in SeqIO.parse(open(fa_filename), 'fasta'):
            zmw = r.id[:r.id.rfind('/')]
            if zmw in zmw_seen:
                continue
            zmw_seen.add(zmw)
            if r.id not in tally:
                unmapped += 1
            else:
                if len(tally[r.id]) == 1:
                    coverages[tally[r.id][0].coverage] += 1
                for gmap_rec in tally[r.id]:
                    for score, seq_exon in zip(gmap_rec.scores, gmap_rec.seq_exons):
                        obs2exp.append((score, score, seq_exon.end-seq_exon.start))
    else:
        for r in SeqIO.parse(open(fq_filename), 'fastq'):
            zmw = r.id[:r.id.rfind('/')]
            if zmw in zmw_seen:
                continue
            zmw_seen.add(zmw)
            if r.id not in tally:
                unmapped += 1
            else:
                for gmap_rec in tally[r.id]:
                    for score, seq_exon in zip(gmap_rec.scores, gmap_rec.seq_exons):
                        obs_acc = np.mean([1-10**-(x/10.) for x in r.letter_annotations['phred_quality'][seq_exon.start:seq_exon.end]])
                        obs2exp.append((obs_acc*100., score, seq_exon.end-seq_exon.start))

    obs2exp = np.array(obs2exp, dtype=[('ObsAccuracy', '>f4'), ('ExpAccuracy', '>f4'), ('Size', '>i4')])
    #os.remove(out_filename)   

    chimera_missed_adapter = 0
    chimera_real = 0
    for v in tally.itervalues():
        if len(v) > 1: # gmap chimera
            if calc_overlap(v[0], v[1]): chimera_missed_adapter += 1
            else: chimera_real += 1

    avg_coverage = sum(k*v for k,v in coverages.iteritems())*1./sum(coverages.itervalues())
    return unmapped, obs2exp, (chimera_missed_adapter,chimera_real,(chimera_missed_adapter+chimera_real)*1./len(tally)), avg_coverage, len(zmw_seen)
        
    
def main(fofn_filename, prefix, random_prob=0.01, num_of_seqs=1000, use_CCS=False, min_seq_len=500):
    out_fa = open(prefix+'.fa', 'w')
    out_fq = open(prefix+'.fq', 'w')
    count = 0
    with open(fofn_filename) as f:
        for line in f:
            print >> sys.stderr, "selecting random sequences from {0}....".format(line.strip())
            count += select_random_sequences(line.strip(), out_fa, out_fq, random_prob, num_of_seqs, use_CCS, min_seq_len)
    out_fa.close()
    out_fq.close()
    
    if count == 0:
        print >> sys.stderr, "Retrieved 0 sequences! :( Abort."
        return None, None
    
    
    print >> sys.stderr, "running gmap on {0}".format(out_fa.name)
    unmapped, obs2exp, chimera_rate, avg_coverage = run_gmap(out_fa.name, out_fq.name)
    
    avg_exp = obs2exp['ExpAccuracy'].mean()
    avg_obs = obs2exp['ObsAccuracy'].mean()
    
    print "Total number of sequences:", count
    print "Total number of unmapped: {0} ({1:.1f}%)".format(unmapped, unmapped*100./count)
    print "Avg. coverage (non-chimeric): {0:.1f}%".format(avg_coverage)
    print "Avg. expected accuracy: {0:.1f}%".format(avg_exp)
    print "Avg. observed accuracy: {0:.1f}%".format(avg_obs)
    print "Chimera rate: {0:.1f}%".format(chimera_rate*100.)

    np.savetxt(prefix+'.obs2exp.txt', obs2exp)
    
    return unmapped, obs2exp
    

def main_eval_only(fa_filename, fq_filename, gmap_filename):
    unmapped, obs2exp, chimera_rate, avg_coverage, zmw_count = run_gmap(fa_filename, fq_filename, gmap_filename)
    
    #count = int(os.popen("grep -c \">\" " + fa_filename).read())
    count = zmw_count

    avg_exp = obs2exp['ExpAccuracy'].mean()
    avg_obs = obs2exp['ObsAccuracy'].mean()
    
    print "Total number of ZMWs:", count
    print "Total number of unmapped: {0} ({1:.1f}%)".format(unmapped, unmapped*100./count)
    print "Avg. coverage (non-chimeric): {0:.1f}%".format(avg_coverage)
    if fq_filename is not None:
        print "Avg. expected accuracy: {0:.1f}%".format(avg_exp)
    print "Avg. observed accuracy: {0:.1f}%".format(avg_obs)
    print "Chimera rate: {0:.1f}%".format(chimera_rate*100.)



def split_gmap_outcome(fa_filename, gmap_filename):
    """
    Split into:
    .gmap_unmapped.fa
    .gmap_chimera.fa
    .gmap_non_chimera.fa
    """
    f_un = open(fa_filename + '.gmap_unmapped.fa', 'w')
    f_is = open(fa_filename + '.gmap_chimera.fa', 'w')
    f_non = open(fa_filename + '.gmap_non_chimera.fa', 'w')

    d = defaultdict(lambda: 0)
    for r in GFF.gmapGFFReader(gmap_filename):
        d[r.seqid] += 1

    for r in SeqIO.parse(open(fa_filename), 'fasta'):
        if r.id not in d: f_un.write(">{0}\n{1}\n".format(r.id, r.seq))
        elif d[r.id] == 1: f_non.write(">{0}\n{1}\n".format(r.id, r.seq))
        else: f_is.write(">{0}\n{1}\n".format(r.id, r.seq))


    f_un.close()
    f_is.close()
    f_non.close()


if __name__ == "__main__":
    global DBNAME
    from argparse import ArgumentParser
    parser = ArgumentParser("Randomly select subreads/CCS reads, run GMAP, plot exp vs obs accuracy")
    parser.add_argument("-i", "--input_fofn", default="input.fofn", help="input fofn")
    parser.add_argument("-p", "--prefix", required=True, help="Output FASTA/FASTQ prefix")
    parser.add_argument("-r", "--random_prob", type=float, default=.01, help="Random selection prob (default: 0.01)")
    parser.add_argument("-n", "--max_seq_per_bash5", type=int, default=1000, help="Max # of seqs per .bas.h5 or .bax.h5 (default: 1000)")
    parser.add_argument("-l", "--min_seq_len", type=int, default=500, help="Minimum seq length (default:500)")
    parser.add_argument("--use_CCS", action='store_true', default=False, help="Use CCS instead of subreads")
    parser.add_argument("--dbname", default='hg19', help="GMAP db name (default: hg19)")

    args = parser.parse_args()
    DBNAME = args.dbname
    main(args.input_fofn, args.prefix, args.random_prob, args.max_seq_per_bash5, args.use_CCS, args.min_seq_len)
