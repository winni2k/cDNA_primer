#!/usr/bin/env python

#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$
__author__ = 'etseng@pacificbiosciences.com'
import os, sys, subprocess
from collections import defaultdict
import numpy as np
from pbcore.io.FastaIO import FastaReader
from pbtools.pbtranscript.io.FastaRandomReader import FastaRandomReader
from argparse import ArgumentParser

__version__ = "Liz-2.1.0"

"""
Precursor to using the smrtpipe pbdagcon
(1) Find the best seed as reference
(2) Align rest to seed
(3) Call pbdagcon
"""

class AlignGraphUtilError(Exception):
    pass

def choose_template_by_blasr(fasta_filename, nproc=8, maxScore=-1000, min_number_reads=1):
    """
    Choose the best template for gcon reference
    Pick the one that has the highest average hit identity to others

    Returns: FastaRecord of selected ref
    """
    fd = FastaRandomReader(fasta_filename)


    out_filename = fasta_filename + '.saln'
    cmd = "blasr -nproc {nproc} -maxScore {score} -maxLCPLength 15 -bestn 10 -nCandidates 50 -m 1 {input} {input} -out {out}".format(\
        nproc=nproc, score=maxScore, input=fasta_filename, out=out_filename)

    if subprocess.check_call(cmd, shell=True) != 0:
        return None

    # blasr -m 1 output format:
    # (0) qName (1) tName (2) qStrand (3) tStrand (4) score (5) percentSimilarity
    # (6) tStart (7) tEnd (8) tLength (9) qStart (10) qEnd (11) qLength (12) nCells
    scores = defaultdict(lambda: [])
    with open(out_filename) as f:
        for line in f:
            raw = line.strip().split()
            qID, tID = raw[0][:raw[0].rfind('/')], raw[1]  # qID gets an extra /0_length
            if qID == tID: continue # self-hit, ignore
            if raw[2] != raw[3]: continue # has to be on same strand
            scores[qID].append(float(raw[5]))  # use identity as the scorer

    # find the one with the highest average alignment similarity
    score_array = []
    for k,v in scores.iteritems():
        score_array.append((np.ceil(np.mean(v)), k))
    if len(score_array) < min_number_reads:
        raise AlignGraphUtilError, ("Not enough number of reads in choose_template_by_blasr {0} < {1}".format(\
            len(score_array), min_number_reads))
    score_array.sort(reverse=True)

    # Liz: find the longest sequence that is within the std deviation of the best score
    best_mean, best_id = score_array[0]
    best_len = len(fd[best_id].sequence)
    for _mean, _id in score_array[1:]:
        if _mean != best_mean: break
        _len = len(fd[_id].sequence)
        if _len > best_len:
            best_id = _id
            best_len = _len

    return fd[best_id]

def make_aln_input_to_ref(fasta_filename, ref_filename, nproc=8):
    """
    Make blasr -m 5 output of input aligned to ref

    Return: alignment filename
    """
    # pbdagcon only takes -m 5 output
    cmd = "blasr {input} {ref} -bestn 1 -nproc {nproc} -m 5 -out {input}.saln.tmp".format(\
        input=fasta_filename, ref=ref_filename, nproc=nproc)

    if subprocess.check_call(cmd, shell=True) != 0:
        raise AlignGraphUtilError, "Unable to align {0} to {1} in make_aln_input_to_ref".format(fasta_filename, ref_filename)

    # trim away anything that is a self-hit or opp strand
    with open(fasta_filename + '.saln', 'w') as f:
        h = open(fasta_filename + '.saln.tmp')
        for line in h:
            raw = line.strip().split()
            # blasr -m 5 output format:
            # (0) qName (1) qLength (2) qStart (3) qEnd (4) qStrand
            # (5) tName (6) tLength (7) tStart (8) tEnd (9) tStrand
            # (10...) score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq
            if raw[0] == raw[5]: continue # self-hit
            if raw[4] != raw[9]: continue # opp strand
            f.write(line)
        h.close()

    os.remove(fasta_filename + '.saln.tmp')
    return f.name


def pbdagcon_wrapper(fasta_filename, output_prefix, consensus_name, nproc=8, maxScore=-1000, min_seq_len=300):
    """
    (1) Find the best seed as reference
    (2) Align rest to seed
    (3) Call pbdagcon
    """
    try:
        ref = choose_template_by_blasr(fasta_filename, nproc=nproc, maxScore=maxScore)

        ref_filename = output_prefix + '_ref.fa'
        with open(ref_filename, 'w') as f:
            f.write(">{0}\n{1}".format(consensus_name, ref.sequence))

        # create alignment file
        aln_filename = make_aln_input_to_ref(fasta_filename, ref_filename, nproc=nproc)

        cons_filename = output_prefix + '.fa'
        # call pbdagcon
        cmd = "pbdagcon -t 0 -m {minlen} -c 1 -j {nproc} {aln} > {out}".format(\
            minlen=min_seq_len, nproc=nproc, aln=aln_filename, out=cons_filename)
        if subprocess.check_call(cmd, shell=True):
            raise AlignGraphUtilError, "Cannot run command:", cmd

    except AlignGraphUtilError:
        # pick the first sequence as reference as a backup plan
        first_seq = FastaReader(fasta_filename).__iter__().next()
        with open(ref_filename, 'w') as f:
            f.write(">{0}_ref\n{1}".format(consensus_name, first_seq.sequence))

def runConsensus(args_list_str):
    parser = ArgumentParser()
    parser.add_argument("input_fasta", help="Input fasta filename")
    parser.add_argument("output_prefix", help="Output filename prefix (ex: g_consensus)")
    parser.add_argument("consensus_id", help="Consensus sequence ID name (ex: consensus)")
    parser.add_argument("--nproc", default=8, type=int, help="Number of processes")
    parser.add_argument("--maxScore", default=-1000, type=int, help="blasr maxScore")
    args_list = args_list_str.split()
    args = parser.parse_args(args_list)
    pbdagcon_wrapper(args.input_fasta, args.output_prefix, args.consensus_id, nproc=args.nproc, maxScore=args.maxScore)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input_fasta", help="Input fasta filename")
    parser.add_argument("output_prefix", help="Output filename prefix (ex: g_consensus)")
    parser.add_argument("consensus_id", help="Consensus sequence ID name (ex: consensus)")
    parser.add_argument("--nproc", default=8, type=int, help="Number of processes")
    parser.add_argument("--maxScore", default=-1000, type=int, help="blasr maxScore")
    args = parser.parse_args()

    pbdagcon_wrapper(args.input_fasta, args.output_prefix, args.consensus_id, nproc=args.nproc, maxScore=args.maxScore)


