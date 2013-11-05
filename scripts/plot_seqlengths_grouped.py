#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

inputs = sys.argv[1].split(",")
names = sys.argv[2].split(",")
styles = sys.argv[3].split(",")
colors = sys.argv[4].split(",")
output = sys.argv[5] + '.seqlengths.png'
range_min = int(sys.argv[6])
range_max = int(sys.argv[7])

fig = plt.figure()
ax1 = fig.add_subplot(111)
bins = (range_max-range_min)/100+1
for input,name,style,color in zip(inputs,names,styles,colors):
    raw = [len(r.seq) for r in SeqIO.parse(open(input), 'fasta')]
    raw = filter(lambda x: range_min<=x<=range_max, raw)

    seqlengths = np.array(raw)
    y,binEdges = np.histogram(seqlengths, bins=bins, normed=True)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

    xnew = np.linspace(bincenters.min(), bincenters.max(), 300)
    ysmooth = spline(bincenters, y, xnew)
    ax1.plot(xnew, ysmooth, style, color=color, label=name)

ax1.legend()
plt.xlabel("Sequence Length")
plt.ylabel("Density")
plt.ylim((0, max(ysmooth)+0.002))
plt.show()
plt.savefig(output,dpi=300)
