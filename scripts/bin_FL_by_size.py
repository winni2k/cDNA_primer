#!/usr/bin/env python
import os
from Bio import SeqIO
from csv import DictReader
seqlen_dict=dict((r.id,len(r.seq)) for r in SeqIO.parse(open('reads_of_insert.fasta'),'fasta'))
bins=[0]*12 # every 500 bp
from collections import defaultdict
bins=defaultdict(lambda: {'fl':0, 'nofl':0})
reader=DictReader(open('reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt'),delimiter='\t')
for r in reader:
    b=seqlen_dict[r['ID'][:r['ID'].rfind('/')]+'/ccs']/500
    if r['5seen']=='1' and r['3seen']=='1' and r['polyAseen']=='1': bins[b]['fl']+=1
    else: bins[b]['nofl']+=1
    
f=open('reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt.bin_by_len.txt','w')
dirname=os.path.abspath('.').split('/')[-1]
f.write(dirname+',')
#for i in xrange(12): print i*500, (i+1)*500-1, bins[i]['fl']*100./(bins[i]['fl']+bins[i]['nofl'])
for i in xrange(12): f.write(str(bins[i]['fl']*100./(bins[i]['fl']+bins[i]['nofl']))+',')
f.write('\n')
