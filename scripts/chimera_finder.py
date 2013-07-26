#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'
import os, sys, multiprocessing, subprocess
from collections import defaultdict
from Bio import SeqIO
from hmmer_wrapper import DOMRecord, sanity_check_phmmer, worker

def sanity_check_primers_for_chimera(primer_filename, p_filename):
    """
    Check that primers are in order of F0, R0, F1, R1, .....
    """
    cur_index = 0
    cur_head = 'F'
    last = None
    f = open(p_filename, 'w')
    for r in SeqIO.parse(open(primer_filename), 'fasta'):
        expected_id = cur_head + str(cur_index)
        if r.id != expected_id:
            print >> sys.stderr, "expected id {0} but got {1}. Bad ID. Abort!".format(expected_id, r.id)
            sys.exit(-1)
        if cur_head == 'R':
            tmp_r = str(last.seq.reverse_complement())
            tmp_seq = str(r.seq)
            if tmp_seq.find(tmp_r) >= 0 or tmp_r.find(tmp_seq) >= 0:
                print >> sys.stderr, "F{0}/R{0} primer pair are reverse-complementarily identical. Adding 'AAAA' in 3' to distinguish".format(cur_index)
                f.write(">{0}\n{1}\n>{2}\nAAAA{3}\n".format(last.id, last.seq, r.id, r.seq))
                f.write(">{0}_revcomp\n{1}\n>{2}_revcomp\n{3}TTTT\n".format(last.id, last.seq.reverse_complement(), r.id, r.seq.reverse_complement()))
            else:
                f.write(">{0}\n{1}\n>{2}\n{3}\n".format(last.id, last.seq, r.id, r.seq))
                f.write(">{0}_revcomp\n{1}\n>{2}_revcomp\n{3}\n".format(last.id, last.seq.reverse_complement(), r.id, r.seq.reverse_complement()))
            cur_index += 1
        last = r
        cur_head = 'R' if cur_head == 'F' else 'F'
    f.close()
    return range(cur_index)

def parse_hmmer_dom_for_chimera(dom_filename, min_score, min_dist_from_end):
    """
    Parses DOM output from phmmer
    Look specifically for hits in the MIDDLE of the sequence as evidence for chimeras
    """
    suspicious_hit = defaultdict(lambda: []) # sid --> list of DOMRecord passing criteria
    with open(dom_filename) as f:
        for line in f:
            if line.startswith('#'): continue # header, ignore
            raw = line.strip().split()
            #pid = raw[0]
            sid = raw[3] # ex: m130212_025538_42161_c100473810150000001823071506131362_s1_p0/45/1687_1987_back
            slen = int(raw[5])
            score = float(raw[13])
            sStart = int(raw[15]) - 1
            sEnd = int(raw[16])
            pStart = int(raw[17]) - 1
            pEnd = int(raw[18])

            # has to be somewhat in the middle
            # and with decent score
            if sStart > min_dist_from_end and sEnd < slen - min_dist_from_end and score >= min_score:
                suspicious_hit[sid].append(DOMRecord(pStart, pEnd, sStart, sEnd, score))
    return suspicious_hit

def remove_chimeras(fasta_filename, suspicious_hits, max_adjacent_hit_distance):
    """
    Output written to <fasta_filename>.non_chimera.fa and <fasta_filename>.is_chimera.fa
    """
    def is_chimera(dom_records):
        """
        Find at least two records such that they are adjacent but NOT entirely overlapping
        """
        if len(dom_records) <= 1: return False
        dom_records.sort(key=lambda x: x.sStart)
        for i in xrange(len(dom_records)-1):
            r1 = dom_records[i]
            s1 = r1.sEnd - r1.sStart
            for j in xrange(i+1, len(dom_records)):
                r2 = dom_records[j]
                s2 = r2.sEnd - r2.sStart
                d = r2.sStart - r1.sEnd
                if d > max_adjacent_hit_distance: break
                elif 0 <= d < max_adjacent_hit_distance:
                    return True
                elif d < 0: # d < 0, has overlap
                    d = abs(d)
                    if d < .5 * s1 and d < .5 * s2: # acceptable overlap for considering them as separate hits (hence is chimera)
                        return True
        return False

    f_chimera = open(fasta_filename + '.is_chimera.fa', 'w')
    f_nonchimera = open(fasta_filename + '.non_chimera.fa', 'w')
    chimera_ids = [sid for (sid, dom_records) in suspicious_hits.iteritems() if is_chimera(dom_records)]
    count_chimera, count_nonchimera = 0, 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        if r.id in chimera_ids:
            f_chimera.write(">{0}\n{1}\n".format(r.id, r.seq))
            count_chimera += 1
        else:
            f_nonchimera.write(">{0}\n{1}\n".format(r.id, r.seq))
            count_nonchimera += 1            
    f_chimera.close()
    f_nonchimera.close()
    print >> sys.stderr, "Number of chimera-to-non-chimera: {0}/{1}".format(count_chimera, count_nonchimera)

def chimera_finder_main(output_dir, primer_filename, fasta_filename, min_dist_from_end=100, max_adjacent_hit_distance=50, cpus=8, min_score=10):
    # find the matrix file PBMATRIX.txt
    matrix_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'PBMATRIX.txt')
    if not os.path.exists(matrix_filename):
        print >> sys.stderr, "Expected matrix file {0} but not found. Abort!".format(matrix_filename)
        sys.exit(-1)
    
    out_filename_hmmer = os.path.join(output_dir, 'hmmer_for_chimera.out')
    if os.path.exists(output_dir):
        if os.path.exists(out_filename_hmmer):
            print >> sys.stderr, "output directory {0} already exists. Running just the primer trimming part.".format(output_dir)
            p_indices = []
            for r in SeqIO.parse(open(os.path.join(output_dir, primer_filename)), 'fasta'):
                if r.id[0] == 'F':
                    p_indices.append(r.id[1:])
        else:
            print >> sys.stderr, "output directory {0} already exists. Abort.".format(output_dir)
            sys.exit(-1)
    else:
        print >> sys.stderr, "checking for phmmer existence..."
        sanity_check_phmmer()
        print >> sys.stderr, "creating output directory {0}....".format(output_dir)
        os.makedirs(output_dir)

        print >> sys.stderr, "checking and copying primer file", primer_filename
        p_filename = os.path.join(output_dir, os.path.basename(primer_filename))
        p_indices = sanity_check_primers_for_chimera(primer_filename, p_filename)

        print >> sys.stderr, "splitting into chunks", fasta_filename
        i = 0
        size = int(os.popen("grep -c \">\" " + fasta_filename).read()) / cpus + 1
        count = 0
        jobs = []
        f_in = open(os.path.join(output_dir, 'in.fa_split'+str(i)), 'w')
        for r in SeqIO.parse(open(fasta_filename), 'fasta'):
            f_in.write(">{0}\n{1}\n".format(r.id, r.seq))
            count += 1
            if count > size:
                f_in.close()
                p = multiprocessing.Process(target=worker, args=(out_filename_hmmer+'_split'+str(i), p_filename, f_in.name, matrix_filename))
                jobs.append((p, out_filename_hmmer+'_split'+str(i)))
                p.start()
                i += 1
                count = 0
                f_in = open(os.path.join(output_dir, 'in.fa_split'+str(i)), 'w')
        f_in.close()
        if count > 0:
            p = multiprocessing.Process(target=worker, args=(out_filename_hmmer+'_split'+str(i), p_filename, f_in.name, matrix_filename))
            jobs.append((p, out_filename_hmmer+'_split'+str(i)))
            p.start()

        for p, out in jobs:
            p.join()
            subprocess.check_call("cat {0} >> {1}".format(out, out_filename_hmmer), shell=True)


    suspicious_hits = parse_hmmer_dom_for_chimera(out_filename_hmmer, min_score, min_dist_from_end)
    remove_chimeras(fasta_filename, suspicious_hits, max_adjacent_hit_distance)
    
    print >> sys.stderr, "Cleaning split files"
    subprocess.check_call("rm -rf {0}/*split*".format(output_dir), shell=True)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="Identify potential chimeras using 5'/3' primers",\
        formatter_class=argparse.RawTextHelpFormatter, \
        description=" This script requires phmmer from HMMER 3.0.\n If the output directory already exists, will skip running phmmer and directory go to primer trimming.\n If you want to re-run HMMER you must first delete the output directory manually.\n Refer to wiki: https://github.com/PacificBiosciences/cDNA_primer/wiki for more details.")

    group1 = parser.add_argument_group("HMMER options")
    group1.add_argument("-p", "--primer_filename", default="primers.fa", help="Primer fasta file")
    group1.add_argument("-i", "--input_filename", default="filtered_subreads.fasta", help="Input fasta file (usually filtered_subreads.fasta or filtered_CCS_subreads.fasta)")
    group1.add_argument("-d", "--directory", default="output", help="Directory to store HMMER output (default: output/)")
    group1.add_argument("--cpus", default=8, type=int, help="Number of CPUs to run HMMER (default: 8)")

    group2 = parser.add_argument_group("Chimera detection options")
    group2.add_argument("--min_dist-from_end", type=int, default=100, help="Minimum distance the primer hit has to be from end of sequence (default: 100)")
    group2.add_argument("--max_adjacent_hit_distance", type=int, default=50, help="Maximum distance between adjacent primer hits to consider as chimera (default: 50)")
    group2.add_argument("--min-score", dest="min_score", type=float, default=10, help="Minimum bit score for primer hit (default: 10)")

    args = parser.parse_args()

    chimera_finder_main(args.directory, args.primer_filename, args.input_filename, args.min_dist_from_end, args.max_adjacent_hit_distance, args.cpus, args.min_score)

    