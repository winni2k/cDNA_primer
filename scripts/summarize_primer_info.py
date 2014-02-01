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
import argparse
from csv import DictReader
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("primer_info_filename", help=".primer_info.txt filename")

args = parser.parse_args()

input = args.primer_info_filename # should be primer_info.txt

# (1) count by subreads
# (2) count by ZMW

num_subread, num_subread_5seen, num_subread_3seen, num_subread_53seen, num_subread_53Aseen = 0, 0, 0, 0, 0
ZMW_5seen = {} # zmw --> list of [True, False, False]....where i-th is i-th subread
ZMW_3seen = {}
ZMW_Aseen = {}
pm_count = defaultdict(lambda: 0)

isCCS = False
with open(input) as f:
    for r in DictReader(f, delimiter='\t'):
        see5 = r['5seen']=='1'
        see3 = r['3seen']=='1'
        seeA = r['polyAseen']=='1'
        if r['primer']!='NA':
            pm_count[r['primer']] += 1
        if r['ID'].count('/') == 1:
            isCCS = True
            zmw = r['ID']
        else:
            zmw = r['ID'][:r['ID'].rfind('/')] # use <movie>/<holeNumber>
        if zmw not in ZMW_5seen:
            ZMW_5seen[zmw] = []
            ZMW_3seen[zmw] = []
            ZMW_Aseen[zmw] = []
        ZMW_5seen[zmw].append(see5)
        ZMW_3seen[zmw].append(see3)
        ZMW_Aseen[zmw].append(seeA)
        num_subread += 1
        num_subread_5seen += see5
        num_subread_3seen += see3
        num_subread_53seen += see5 and see3
        num_subread_53Aseen += see5 and see3 and seeA

print "------ 5' primer seen summary ---- "
if not isCCS:
    print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_5seen, num_subread, num_subread_5seen*100./num_subread)
tmp = sum(any(x) for x in ZMW_5seen.itervalues())
num_ZMW = len(ZMW_5seen)
assert num_ZMW == len(ZMW_3seen)
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum((len(x)>0 and x[0]) for x in ZMW_5seen.itervalues())
if not isCCS:
    print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)

print "------ 3' primer seen summary ---- "
if not isCCS:
    print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_3seen, num_subread, num_subread_3seen*100./num_subread)
tmp = sum(any(x) for x in ZMW_3seen.itervalues())
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum((len(x)>0 and x[0]) for x in ZMW_3seen.itervalues())
if not isCCS:
    print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)

print "------ 5'&3' primer seen summary ---- "
if not isCCS:
    print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_53seen, num_subread, num_subread_53seen*100./num_subread)
tmp = 0
for zmw,x in  ZMW_5seen.iteritems():
    for i in xrange(len(x)):
        if x[i] and ZMW_3seen[zmw][i]:
            tmp += 1
            break
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum([(len(x)>0 and x[0] and len(ZMW_3seen[zmw])>0 and ZMW_3seen[zmw][0]) for (zmw,x) in ZMW_5seen.iteritems()])
if not isCCS:
    print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)


print "------ 5'&3'&polyA primer seen summary ---- "
if not isCCS:
    print "Per subread: {0}/{1} ({2:.1f}%)".format(num_subread_53Aseen, num_subread, num_subread_53Aseen*100./num_subread)
tmp = 0
for zmw,x in  ZMW_5seen.iteritems():
    for i in xrange(len(x)):
        if x[i] and ZMW_3seen[zmw][i] and ZMW_Aseen[zmw][i]:
            tmp += 1
            break
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum([(len(x)>0 and x[0] and len(ZMW_3seen[zmw])>0 and ZMW_3seen[zmw][0] and len(ZMW_Aseen[zmw])>0 and ZMW_Aseen[zmw][0]) for (zmw,x) in ZMW_5seen.iteritems()])
if not isCCS:
    print "Per ZMW first-pass: {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)


print "------ Primer Match breakdown ---- "
keys = pm_count.keys()
keys.sort()
total = sum(pm_count.itervalues())
for k in keys:
    print "F{0}/R{0}: {1} ({2:.1f}%)".format(k, pm_count[k], pm_count[k]*100./total)
