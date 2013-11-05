#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

inputs = sys.argv[1].split(",")
names = sys.argv[2].split(",")
styles = sys.argv[3].split(",")
colors = sys.argv[4].split(",")
output = sys.argv[5] + '.FL_by_bin.png'

fig = plt.figure()
ax1 = fig.add_subplot(111)
for input,name,style,color in zip(inputs,names,styles,colors):
    y = map(float, open(input).readline().strip().split(',')[1:-1])
    x = xrange(0, 6000, 500)
    ax1.plot(x, y, style, color=color, label=name)

ax1.legend()
plt.xlabel("Sequence Length Range")
plt.ylabel("FL %")
plt.ylim((0, 100))
plt.show()
plt.savefig(output,dpi=300)
