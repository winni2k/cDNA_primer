import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import BioReaders

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
    """
    Return dict of refID --> count of ZMW hits
    """
    for r in hit_by_zmw.itervalues():
        if r.sCoverage >= min_sCov and r.qCoverage >= min_qCov:
            ref_tally[r.sID] += 1
            

def draw_2dhist(hit, output_prefix, feat_func=lambda x: (x.sLen, x.qLen), filter_func=lambda x: True, xlab='Reference Length', ylab='Query length', title=''):
    """
    Heatmap plot of query (CCS/non-CCS) seq length VS reference length (or whatever is given)
    for which there is alignment
    (NOTE: this is drawing the length of the input seq, not the aligned length!)
    """
    x = []
    y = []
    for r in hit.itervalues():
        if filter_func(r):
            a, b = feat_func(r) 
            x.append(a)
            y.append(b)
    x = np.array(x)
    y = np.array(y)
    
    fig = plt.figure(dpi=300, figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.hexbin(x, y, bins='log', edgecolors='none', cmap=matplotlib.cm.jet)
    
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.savefig(output_prefix+'.png', format='png')


def draw_coverage(hit, output_prefix, filter_func=lambda x: x.qCoverage>=.9 and 1000<=x.sLen<=2000, title="5'-3' reference coverage (qCov>=90%), ref length 1000-2000bp"):
    """
    5'-3' reference coverage, normalized by length
    """
    data = np.zeros(101) # divide by 0.01
    for r in hit.itervalues():
        if filter_func(r):
            s = r.sStart * 100 / r.sLen
            e = r.sEnd * 100 / r.sLen
            data[s:(e+1)] += 1
    
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.fill_between(np.arange(0,1.01,0.01), data, facecolor='#3399FF', alpha=.7)
    ax.set_xlabel("Reference length normalized")
    ax.set_ylabel("Count")
    ax.set_title(title)
    fig.savefig(output_prefix+'.png', format='png')
    
    

    
