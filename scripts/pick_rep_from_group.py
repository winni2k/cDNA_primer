#!/usr/bin/bash
import os, sys
from SeqReaders import LazyFastaReader
from pbcore.io.FastaIO import FastaWriter

def main(group_file, fasta_file1, fasta_file2):
    print >> sys.stderr, "Lazy indexing of file {0}...".format(fasta_file1)
    fd1 = LazyFastaReader(fasta_file1)
    print >> sys.stderr, "Lazy indexing of file {0}...".format(fasta_file2)
    fd2 = LazyFastaReader(fasta_file2)
    print >> sys.stderr, "Reading group file {0}...".format(group_file)
    fout1 = FastaWriter(fasta_file1 + '.rep.fa')
    fout2 = FastaWriter(fasta_file2 + '.rep.fa')
    with open(group_file) as f:
        for line in f:
            pid, members = line.strip().split('\t')
            picked = members.split(',')[0]
            print >> sys.stderr, "pick {0} for {1}".format(picked, pid)
            fout1.writeRecord("{0}|{1}/1".format(pid, picked), fd1[picked+'/1'].sequence)
            fout2.writeRecord("{0}|{1}/2".format(pid, picked), fd2[picked+'/2'].sequence)
    fout1.close()
    fout2.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
