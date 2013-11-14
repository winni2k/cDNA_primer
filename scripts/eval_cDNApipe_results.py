#!/usr/bin/env python
import os, sys
import numpy as np
from collections import defaultdict
from Bio import SeqIO

import eval_accuracy_by_gmap
import eval_refmap_by_blasr

def get_zmw(filename):
    """
    Returns dict of {zmw} --> count of sequences (for RoI, must be exactly 1)
    """
    d = defaultdict(lambda: 0)
    for r in SeqIO.parse(open(filename), 'fasta'):
        if r.id.endswith('/ccs'): #  old version of RoI id
            zmw = r.id[:-4]
        elif r.id.count('/') == 1: # new version of RoI
            zmw = r.id
        elif r.id.count('/') == 2: # subread
            zmw = r.id[:r.id.rfind('/')]
        else:
            raise Exception, "Unrecognized ID format: {0}".format(r.id)
        d[zmw] += 1
    return d

def get_readlen(filename): return dict((r.id, len(r.seq)) for r in SeqIO.parse(open(filename), 'fasta'))

def summarize_RoI(roi_filename='reads_of_insert.fasta'):
    """
    Output:
    1. Number of RoI
    2. Mean RoI readlength
    """
    d_roi = get_zmw(roi_filename)
    r_roi = get_readlen(roi_filename)
    
    a = len(d_roi)
    print("Number of RoI reads: {0}".format(len(d_roi)))
    print("Avg. RoI readlength: {0}".format(sum(r_roi.itervalues())/a))
    
def summarize_chimera(roi_prefix='reads_of_insert.53Aseen_trimmed_changeid.fa'):
    """
    Output: % of artificial chimeras in RoI
    """
    d_roi_is = get_zmw(roi_prefix + '.is_chimera.fa')
    d_roi_not = get_zmw(roi_prefix + '.non_chimera.fa')
    
    a = len(d_roi_is)
    b = a + len(d_roi_not)
    print("% of artificial chimeras in RoI: {0}/{1} ({2:.1f}%)").format(a, b, a*100./b)
    
def summarize_gmap(roi_prefix='reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa'):
    roi_unmapped, roi_obs2exp, (roi_chimera_missed,roi_chimera_real,roi_chimera_rate), roi_avg_coverage, roi_zmw_count = eval_accuracy_by_gmap.run_gmap(roi_prefix, None, roi_prefix+'.gff')
        
    roi_avg_obs = roi_obs2exp['ObsAccuracy'].mean()
    
    print "----- GMAP ------"
    print "Total number of ZMWs:", roi_zmw_count
    print "Total number of unmapped: {0} ({1:.1f}%)".format(roi_unmapped, roi_unmapped*100./roi_zmw_count)
    print "Avg. coverage: {0:.1f}%".format(roi_avg_coverage)
    print "Avg. observed accuracy: {0:.1f}%".format(roi_avg_obs)
    print "Chimera rate: {0:.1f}%".format(roi_chimera_rate*100.)
    print "Chimera from overpriming: {0:.1f}%".format(roi_chimera_missed*100./(roi_chimera_missed+roi_chimera_real))

    eval_accuracy_by_gmap.split_gmap_outcome(roi_prefix, roi_prefix+'.gff')

def parse_blasr(sam_filename, ref_fasta_filename):
    """
    Return dict of ZMW --> best r by maximizing sCov
    """
    hit = {}
    for r in BioReaders.BLASRSAMReader(sam_filename, True, ref_fasta_filename):
        zmw = r.qID[:r.qID.find('/', r.qID.find('/')+1)]
        if zmw not in hit or r.sCoverage >= hit[zmw].sCoverage:
            hit[zmw] = r
    return hit

def tally_ref_hits(hit_by_zmw, ref_tally, min_sCov, min_qCov):
    for r in hit_by_zmw.itervalues():
        if r.sCoverage >= min_sCov and r.qCoverage >= min_qCov:
            ref_tally[r.sID] += 1

def summarize_blasr(roi_prefix='reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa', ref_fasta_filename='/home/UNIXHOME/etseng/share/gencode/gencode.v15.pc_transcripts.non_redundant_good.fa.nonredundant.fasta', ref_size=None, png_name_prefix=''):
    """
    1. Number of aligned ZMWs / refs
    2. Number of well-aligned (90% qCov, sCov) ZMWs / refs
    3. Abundance distribution of well-aligned refs
    """
    ref_len_dict = dict((r.id,len(r.seq)) for r in SeqIO.parse(open(ref_fasta_filename),'fasta'))
    hit_roi = eval_refmap_by_blasr.parse_blasr(roi_prefix+'.blasr.sam', ref_len_dict)

    tally0 = defaultdict(lambda: 0)
    tally90 = defaultdict(lambda: 0)

    eval_refmap_by_blasr.tally_ref_hits(hit_roi, tally0, 0., 0.)
    eval_refmap_by_blasr.tally_ref_hits(hit_roi, tally90, 0.9, 0.9)

    ref_count = int(os.popen("grep -c \">\" " + ref_fasta_filename).read().strip())
    zmw_count = len(get_zmw(roi_prefix)) 

    a = len(tally0)
    b = ref_count
    c = len(hit_roi)
    d = zmw_count
    print "----- BLASR ------"
    print("Number of aligned refs: {0}/{1} ({2:.1f}%)".format(a, b, 100.*a/b))
    print("Number of aligned ZMWs: {0}/{1} ({2:.1f}%)".format(c, d, 100.*c/d))
    a = len(tally90)
    c = sum(r.sCoverage>=.9 and r.qCoverage>=.9 for r in hit_roi.itervalues())
    print("Number of well-aligned refs: {0}/{1} ({2:.1f}%)".format(a, b, 100.*a/b))
    print("Number of well-aligned ZMWs: {0}/{1} ({2:.1f}%)".format(c, d, 100.*c/d))

    # draw plots
    hit = hit_roi

    eval_refmap_by_blasr.draw_2dhist(hit, png_name_prefix+'.roi.sLen_vs_qLen', feat_func=lambda x: (x.sLen, x.qLen), filter_func=lambda x: True, xlab='Reference Length', ylab='Query length')
    eval_refmap_by_blasr.draw_2dhist(hit, png_name_prefix+'.roi.sLen_vs_qLen_well', feat_func=lambda x: (x.sLen, x.qLen), filter_func=lambda x: x.qCoverage>=.9 and x.sCoverage>=.9, xlab='Reference Length', ylab='Query length')
    eval_refmap_by_blasr.draw_2dhist(hit, png_name_prefix+'.roi.sCov_vs_qCov', feat_func=lambda x: (x.sCoverage, x.qCoverage), filter_func=lambda x: True, xlab='Reference Coverage', ylab='Query Coverage')
    eval_refmap_by_blasr.draw_2dhist(hit, png_name_prefix+'.roi.sLen_vs_sCov', feat_func=lambda x: (x.sLen, x.sCoverage), filter_func=lambda x: True, xlab='Reference Length', ylab='Reference Coverage')

    if ref_size is None:
        eval_refmap_by_blasr.draw_coverage(hit, png_name_prefix+'.roi.ref_coverage', filter_func=lambda x: x.qCoverage>=.9, title="5'-3' reference coverage (qCov>=90%)")
    else:
        a, b = ref_size
        eval_refmap_by_blasr.draw_coverage(hit, png_name_prefix+'.roi.ref_coverage', filter_func=lambda x: x.qCoverage>=.9 and a<=x.sLen<=b, title="5'-3' reference coverage (qCov>=90%), ref length {0}-{1} bp".format(a,b))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--ref_fasta_filename", default='/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/Gencode15/sequence/Gencode15.fasta')
    parser.add_argument("--ref_size", default=None)
    parser.add_argument("--skip-GMAP", default=False, dest="skip_GMAP", action="store_true", help="Skip GMAP evaluation")
    parser.add_argument("--skip-BLASR", default=False, dest="skip_BLASR", action="store_true", help="Skip BLASR evaluation")
    
    args = parser.parse_args()


    name = os.path.basename(os.getcwd())

    summarize_RoI()
    summarize_chimera()
    if not args.skip_GMAP:
        summarize_gmap()
    if not args.skip_BLASR:
        summarize_blasr(ref_fasta_filename=args.ref_fasta_filename, ref_size=args.ref_size, png_name_prefix=name)
    
    # plot FL seq length
    cmd = "plot_seqlengths_grouped.py reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa {name} - blue {name}.FL_seqlength 0 6000".format(name=name)
    print >> sys.stderr, "CMD:", cmd
    assert os.system(cmd) == 0
    # plot FL % bin by length
    cmd = "bin_FL_by_size.py"
    print >> sys.stderr, "CMD:", cmd
    assert os.system(cmd) == 0
    cmd = "plot_FL_by_in.py reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt.bin_by_len.txt {name} - purple {name}".format(name=name)
    print >> sys.stderr, "CMD:", cmd
    assert os.system(cmd) == 0
