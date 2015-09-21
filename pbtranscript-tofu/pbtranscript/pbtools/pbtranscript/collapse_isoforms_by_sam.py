#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'
#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$
import os, sys
from pbtools.pbtranscript.Utils import check_ids_unique
from pbtools.pbtranscript.io.SeqReaders import LazyFastaReader, LazyFastqReader
from pbtools.pbtranscript.branch import branch_simple2
from pbcore.io.FastaIO import FastaWriter
from pbcore.io.FastqIO import FastqWriter
from pbtools.pbtranscript.counting import compare_junctions
from pbtools.pbtranscript.io import GFF

def pick_rep(fa_fq_filename, gff_filename, group_filename, output_filename, is_fq=False, pick_least_err_instead=False, bad_gff_filename=None):
    """
    For each group, select the representative record

    If is FASTA file (is_fa False) -- then always pick the longest one
    If is FASTQ file (is_fq True) -- then 
          If pick_least_err_instead is True, pick the one w/ least number of expected base errors
          Else, pick the longest one
    """
    if is_fq:
        fd = LazyFastqReader(fa_fq_filename)
        fout = FastqWriter(output_filename)
    else:
        fd = LazyFastaReader(fa_fq_filename)
        fout = FastaWriter(output_filename)

    coords = {}
    for line in open(gff_filename):
        # ex: chr1    PacBio  transcript      27567   29336   .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";
        raw = line.strip().split('\t')
        if raw[2] == 'transcript': 
            tid = raw[-1].split('; ')[1].split()[1][1:-2]
            coords[tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

    if bad_gff_filename is not None:
        for line in open(bad_gff_filename):
            raw = line.strip().split('\t')
            if raw[2] == 'transcript':
                tid = raw[-1].split('; ')[1].split()[1][1:-2]
                coords[tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

    for line in open(group_filename):
        pb_id, members = line.strip().split('\t')
        print >> sys.stderr, "Picking representative sequence for", pb_id
        best_id = None
        best_seq = None
        best_qual = None
        best_err = 9999999
        err = 9999999
        max_len = 0
        for x in members.split(','):
            if is_fq and pick_least_err_instead:
                err = sum(i**-(i/10.) for i in fd[x].quality)
            if (is_fq and pick_least_err_instead and err < best_err) or ((not is_fq or not pick_least_err_instead) and len(fd[x].sequence) >= max_len):
                best_id = x
                best_seq = fd[x].sequence
                if is_fq:
                    best_qual = fd[x].quality
                    best_err = err
                max_len = len(fd[x].sequence)

        _id_ = "{0}|{1}|{2}".format(pb_id, coords[pb_id], best_id)
        _seq_ = best_seq
        if is_fq:
            fout.writeRecord(_id_, _seq_, best_qual)
        else:
            fout.writeRecord(_id_, _seq_)
    fout.close()

from collections import defaultdict
from bx.intervals import IntervalTree

def collapse_fuzzy_junctions(gff_filename, group_filename, allow_extra_5exon, internal_fuzzy_max_dist):
    def get_fl_from_id(members):
        # ex: 13cycle_1Mag1Diff|i0HQ_SIRV_1d1m|c139597/f1p0/178
        return sum(int(_id.split('/')[1].split('p')[0][1:]) for _id in members)

    def can_merge(m, r1, r2):
        if m == 'exact':
            return True
        else:
            if not allow_extra_5exon:
                return False
        # below is continued only if (a) is 'subset' or 'super' AND (b) allow_extra_5exon is True
        if m == 'subset':
            r1, r2 = r2, r1 #  rotate so r1 is always the longer one
        if m == 'super' or m == 'subset':
            n2 = len(r2.ref_exons)
            # check that (a) r1 and r2 end on same 3' exon, that is the last acceptor site agrees
            # AND (b) the 5' start of r2 is sandwiched between the matching r1 exon coordinates
            if r1.strand == '+':
                return abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[-n2].start <= r2.ref_exons[0].start < r1.ref_exons[-n2].end
            else:
                return abs(r1.ref_exons[0].end - r2.ref_exons[0].end) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[n2-1].start <= r2.ref_exons[-1].end < r1.ref_exons[n2].end
        return False

    d = {}
    recs = defaultdict(lambda: {'+':IntervalTree(), '-':IntervalTree()}) # chr --> strand --> tree
    fuzzy_match = defaultdict(lambda: [])
    for r in GFF.collapseGFFReader(gff_filename):
        d[r.seqid] = r
        has_match = False
        r.segments = r.ref_exons
        for r2 in recs[r.chr][r.strand].find(r.start, r.end):
            r2.segments = r2.ref_exons
            m = compare_junctions.compare_junctions(r, r2, internal_fuzzy_max_dist=internal_fuzzy_max_dist)
            if can_merge(m, r, r2):
                fuzzy_match[r2.seqid].append(r.seqid)
                has_match = True
                break
        if not has_match:
            recs[r.chr][r.strand].insert(r.start, r.end, r)
            fuzzy_match[r.seqid] = [r.seqid]

    group_info = {}
    with open(group_filename) as f:
        for line in f:
            pbid, members = line.strip().split('\t')
            group_info[pbid] = [x for x in members.split(',')]

    # pick for each fuzzy group the one that has the most exons (if tie, then most FL)
    keys = fuzzy_match.keys()
    keys.sort(key=lambda x: map(int, x.split('.')[1:]))
    f_gff = open(gff_filename+'.fuzzy', 'w')
    f_group = open(group_filename+'.fuzzy', 'w')
    for k in keys:
        all_members = []
        best_pbid, best_size, best_num_exons = fuzzy_match[k][0], len(group_info[fuzzy_match[k][0]]), len(d[fuzzy_match[k][0]].ref_exons)
        all_members += group_info[fuzzy_match[k][0]]
        for pbid in fuzzy_match[k][1:]:
            _size = get_fl_from_id(group_info[pbid])
            _num_exons = len(d[pbid].ref_exons)
            all_members += group_info[pbid]
            if _num_exons > best_num_exons or (_num_exons == best_num_exons and _size > best_size):
                best_pbid, best_size, best_num_exons = pbid, _size, _num_exons
        GFF.write_collapseGFF_format(f_gff, d[best_pbid])
        f_group.write("{0}\t{1}\n".format(best_pbid, ",".join(all_members)))
    f_gff.close()
    f_group.close()

    return fuzzy_match

    
def main(args):    
    if not os.path.exists(args.input):
        print >> sys.stderr, "Input file {0} does not exist. Abort.".format(args.fasta)
        sys.exit(-1)
    
    if not os.path.exists(args.sam):
        print >> sys.stderr, "SAM file {0} does not exist. Abort.".format(args.sam)
        sys.exit(-1)

    # check for duplicate IDs
    check_ids_unique(args.input, is_fq=args.fq)
    
    ignored_fout = open(args.prefix + '.ignored_ids.txt', 'w')

    if args.flnc_coverage > 0:
        f_good = open(args.prefix + '.collapsed.good.gff', 'w')
        f_bad = open(args.prefix + '.collapsed.bad.gff', 'w')
        cov_threshold = args.flnc_coverage
    else:
        f_good = open(args.prefix + '.collapsed.gff', 'w')
        f_bad = f_good
        cov_threshold = 1
    f_txt = open(args.prefix + '.collapsed.group.txt', 'w')
    
    b = branch_simple2.BranchSimple(args.input, cov_threshold=cov_threshold, min_aln_coverage=args.min_aln_coverage, min_aln_identity=args.min_aln_identity, is_fq=args.fq)
    iter = b.iter_gmap_sam(args.sam, ignored_fout)
    for recs in iter:
        for v in recs.itervalues():
            if len(v) > 0: b.process_records(v, args.allow_extra_5exon, False, f_good, f_bad, f_txt)
    
    ignored_fout.close()
    f_good.close()
    try:
        f_bad.close()
    except:
        pass
    f_txt.close()

    if args.max_fuzzy_junction > 0: # need to further collapse those that have fuzzy junctions!
        collapse_fuzzy_junctions(f_good.name, f_txt.name, args.allow_extra_5exon, internal_fuzzy_max_dist=args.max_fuzzy_junction)
        os.rename(f_good.name, f_good.name+'.unfuzzy')
        os.rename(f_txt.name, f_txt.name+'.unfuzzy')
        os.rename(f_good.name+'.fuzzy', f_good.name)
        os.rename(f_txt.name+'.fuzzy', f_txt.name)

    if args.fq:
        outfile = args.prefix+".collapsed.rep.fq"
    else:
        outfile = args.prefix+".collapsed.rep.fa"
    if args.allow_extra_5exon: # 5merge, pick longest
        pick_rep(args.input, f_good.name, f_txt.name, outfile, is_fq=args.fq, pick_least_err_instead=False, bad_gff_filename=f_bad.name)
    else:
        pick_rep(args.input, f_good.name, f_txt.name, outfile, is_fq=args.fq, pick_least_err_instead=True, bad_gff_filename=f_bad.name)
    
    print >> sys.stderr, "Ignored IDs written to:", ignored_fout.name
    print >> sys.stderr, "Output written to:"
    print >> sys.stderr, f_good.name
    print >> sys.stderr, f_txt.name
    print >> sys.stderr, outfile
    print >> sys.stderr, args


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--input", help="Input FA/FQ filename")
    parser.add_argument("--fq", default=False, action="store_true", help="Input is a fastq file (default is fasta)")
    parser.add_argument("-s", "--sam", required=True, help="Sorted GMAP SAM filename")
    parser.add_argument("-o", "--prefix", required=True, help="Output filename prefix")
    parser.add_argument("-c", "--min-coverage", dest="min_aln_coverage", type=float, default=.99, help="Minimum alignment coverage (default: 0.99)")
    parser.add_argument("-i", "--min-identity", dest="min_aln_identity", type=float, default=.95, help="Minimum alignment identity (default: 0.95)")
    parser.add_argument("--max_fuzzy_junction", default=5, type=int, help="Max fuzzy junction dist (default: 5 bp)")
    parser.add_argument("--flnc_coverage", dest="flnc_coverage", type=int, default=-1, help="Minimum # of FLNC reads, only use this for aligned FLNC reads, otherwise results undefined!")
    parser.add_argument("--dun-merge-5-shorter", action="store_false", dest="allow_extra_5exon", default=True, help="Don't collapse shorter 5' transcripts (default: turned off)")
    
    args = parser.parse_args()
    
    main(args)

