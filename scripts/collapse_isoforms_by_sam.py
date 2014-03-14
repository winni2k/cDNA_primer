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
import branch_simple2
from SeqReaders import LazyFastaReader
from pbcore.io.FastaIO import FastaWriter

def pick_longest_rep(fasta_filename, gff_filename, group_filename, output_filename):
    """
    For each group, select the representative record to be the longest
    """
    fastad = LazyFastaReader(fasta_filename)
    fout = FastaWriter(output_filename)

    coords = {}
    for line in open(gff_filename):
        # ex: chr1    PacBio  transcript      27567   29336   .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";
        raw = line.strip().split('\t')
        if raw[2] == 'transcript': 
            tid = raw[-1].split('; ')[1].split()[1][1:-2]
            coords[tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

    for line in open(group_filename):
        pb_id, members = line.strip().split('\t')
        best_id = None
        best_seq = None
        max_len = 0
        for x in members.split(','):
            if len(fastad[x].sequence) >= max_len:
                best_id = x
                best_seq = fastad[x].sequence
                max_len = len(fastad[x].sequence)
        fout.writeRecord("{0}|{1}|{2}".format(pb_id, coords[pb_id], best_id), best_seq)
    fout.close()
    
def main(args):    
    if not os.path.exists(args.fasta):
        print >> sys.stderr, "Fasta file {0} does not exist. Abort.".format(args.fasta)
        sys.exit(-1)
    
    if not os.path.exists(args.sam):
        print >> sys.stderr, "SAM file {0} does not exist. Abort.".format(args.sam)
        sys.exit(-1)
    
    ignored_fout = open(args.prefix + '.ignored_ids.txt', 'w')
    f_gff = open(args.prefix + '.collapsed.gff', 'w')
    f_txt = open(args.prefix + '.collapsed.group.txt', 'w')
    
    b = branch_simple2.BranchSimple(args.fasta, cov_threshold=1, min_aln_coverage=args.min_aln_coverage, min_aln_identity=args.min_aln_identity)
    iter = b.iter_gmap_sam(args.sam, ignored_fout)
    for recs in iter:
        for v in recs.itervalues():
            if len(v) > 0: b.process_records(v, args.allow_extra_5exon, False, f_gff, f_gff, f_txt)
    
    ignored_fout.close()
    f_gff.close()
    f_txt.close()
    
    pick_longest_rep(args.fasta, f_gff.name, f_txt.name, args.prefix+".collapsed.longest_rep.fa")
    
    print >> sys.stderr, "Ignored IDs written to:", ignored_fout.name
    print >> sys.stderr, "Output written to:"
    print >> sys.stderr, f_gff.name
    print >> sys.stderr, f_txt.name
    print >> sys.stderr, args.prefix+".collapsed.longest_rep.fa"



if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True, help="Fasta filename")
    parser.add_argument("-s", "--sam", required=True, help="Sorted GMAP SAM filename")
    parser.add_argument("-o", "--prefix", required=True, help="Output filename prefix")
    parser.add_argument("-c", "--min-coverage", dest="min_aln_coverage", default=.99, help="Minimum alignment coverage (default: 0.99)")
    parser.add_argument("-i", "--min-identity", dest="min_aln_identity", default=.85, help="Minimum alignment identity (default: 0.85)")
    parser.add_argument("--dun-merge-5-shorter", action="store_false", dest="allow_extra_5exon", default=True, help="Don't collapse shorter 5' transcripts (default: turned off)")
    
    args = parser.parse_args()
    
    main(args)

