#!/usr/bin/env python
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

num_read, num_read_fiveseen, num_read_threeseen, num_read_5threeseen, num_read_53Aseen = 0, 0, 0, 0, 0
ZMW_fiveseen = {} # zmw --> list of [True, False, False]....where i-th is i-th subread
ZMW_threeseen = {}
ZMW_Aseen = {}
pm_count = defaultdict(lambda: 0)

with open(input) as f:
    for r in DictReader(f, delimiter='\t'):
        see5 = r['fiveseen']=='1'
        see3 = r['threeseen']=='1'
        seeA = r['polyAseen']=='1'
        if r['primer']!='NA':
            pm_count[r['primer']] += 1
        if r['ID'].count('/') == 1:
            isCCS = True
            zmw = r['ID']
        else:
            zmw = r['ID'][:r['ID'].rfind('/')] # use <movie>/<holeNumber>
        if zmw not in ZMW_fiveseen:
            ZMW_fiveseen[zmw] = []
            ZMW_threeseen[zmw] = []
            ZMW_Aseen[zmw] = []
        ZMW_fiveseen[zmw].append(see5)
        ZMW_threeseen[zmw].append(see3)
        ZMW_Aseen[zmw].append(seeA)
        num_read += 1
        num_read_fiveseen += see5
        num_read_threeseen += see3
        num_read_5threeseen += see5 and see3
        num_read_53Aseen += see5 and see3 and seeA

print "------ 5' primer seen summary ---- "
tmp = sum(any(x) for x in ZMW_fiveseen.itervalues())
num_ZMW = len(ZMW_fiveseen)
assert num_ZMW == len(ZMW_threeseen)
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum((len(x)>0 and x[0]) for x in ZMW_fiveseen.itervalues())

print "------ 3' primer seen summary ---- "
tmp = sum(any(x) for x in ZMW_threeseen.itervalues())
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum((len(x)>0 and x[0]) for x in ZMW_threeseen.itervalues())

print "------ 5'&3' primer seen summary ---- "
tmp = 0
for zmw,x in  ZMW_fiveseen.iteritems():
    for i in xrange(len(x)):
        if x[i] and ZMW_threeseen[zmw][i]:
            tmp += 1
            break
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum([(len(x)>0 and x[0] and len(ZMW_threeseen[zmw])>0 and ZMW_threeseen[zmw][0]) for (zmw,x) in ZMW_fiveseen.iteritems()])

print "------ 5'&3'&polyA primer seen summary ---- "
tmp = 0
for zmw,x in  ZMW_fiveseen.iteritems():
    for i in xrange(len(x)):
        if x[i] and ZMW_threeseen[zmw][i] and ZMW_Aseen[zmw][i]:
            tmp += 1
            break
print "Per ZMW:     {0}/{1} ({2:.1f}%)".format(tmp, num_ZMW, tmp*100./num_ZMW)
tmp = sum([(len(x)>0 and x[0] and len(ZMW_threeseen[zmw])>0 and ZMW_threeseen[zmw][0] and len(ZMW_Aseen[zmw])>0 and ZMW_Aseen[zmw][0]) for (zmw,x) in ZMW_fiveseen.iteritems()])

print "------ Primer Match breakdown ---- "
keys = pm_count.keys()
keys.sort()
total = sum(pm_count.itervalues())
for k in keys:
    print "F{0}/R{0}: {1} ({2:.1f}%)".format(k, pm_count[k], pm_count[k]*100./total)
