#!/usr/bin/env python
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
from collections import defaultdict
from csv import DictReader
from Bio import SeqIO
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-d", dest="hmmer_dir", default="output", help="hmmer output directory from running hmmer_wrapper.py, [default: output]")
parser.add_argument("-i", dest="input_fasta", default="reads_of_insert.fasta", help="Input file that was used to run hmmer_wrapper.py, [default: reads_of_insert.fasta]")
parser.add_argument("-o", dest="output_fasta", default="isoseq_nfl.fasta", help="Ouptut filename (default: isoseq_nfl.fasta)")
parser.add_argument("--min-seqlen", dest="min_seqlen", type=int, default=100, help="Minimum sequence length (default: 100)")
parser.add_argument("--fl-has-no-polyA", dest="ignore_polyA", action="store_true", default=False, help="Use this option if your FL reads have no polyA tail")

args = parser.parse_args()

if not os.path.exists(args.hmmer_dir):
    print >> sys.stderr, "Directory {0} does not exist! Please indicate the directory in which the hmmer_wrapper.py output directory was stored, it usually contains a file called hmmer.out".format(args.hmmer_dir)
    sys.exit(-1)
if not os.path.exists(args.input_fasta):
    print >> sys.stderr, "Fasta file {0} does not exist. Abort!".format(args.input_fasta)
    sys.exit(-1)

cmd = "hmmer_wrapper.py -d {d} -i {i} --output-anyway --change-seqid -p primers.front_end.fa --hmmer-out hmmer.front_end.dom -o reads_of_insert.output_anyway_changeid.fa --min-seqlen {m}".format(d=args.hmmer_dir, i=args.input_fasta, m=args.min_seqlen)
if subprocess.check_call(cmd, shell=True)!=0:
    print >> sys.stderr, "TROUBLE RUNNING COMMAND:", cmd
    sys.exit(-1)

reader=DictReader(open('reads_of_insert.output_anyway_changeid.fa.primer_info.txt'),delimiter='\t')
d=defaultdict(lambda: [])
for r in reader:
    d[r['ID'][:r['ID'].rfind('/')]].append(r)
    
reader=SeqIO.to_dict(SeqIO.parse(open('reads_of_insert.output_anyway_changeid.fa'), 'fasta'))
g=(lambda lst: any((x['fiveseen']=='1' and x['threeseen']=='1' and (args.ignore_polyA or x['polyAseen']=='1')) for x in lst))
f = open(args.output_fasta, 'w')
for v in d.itervalues():
    if g(v): continue
    for x in v:
        try:
            r = reader[x['ID']]
            f.write(">{0}\n{1}\n".format(r.id,r.seq))
        except:
            print >> sys.stderr, "ignore", x['ID']
        
f.close()

os.rename('reads_of_insert.output_anyway_changeid.fa.primer_info.txt', args.output_fasta+'.primer_info.txt')
os.remove('reads_of_insert.output_anyway_changeid.fa') # remove this since unlikely I'll use it
