#!/usr/bin/env python
import os, sys

def check_file_or_run(filename):
    def func(f):
        def g(bash_f, args):
            if args.force or not os.path.exists(filename):
                f(bash_f, args)
            else:
                print >> sys.stderr, "{0} already exists. No need to run {1}.".format(filename, f.__name__)
        return g
    return func

@check_file_or_run('reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt.summary')
def get_primer_info(bash_f, args):
    """
    params: --change-seqid, --must-see-polyA --min-seqlen 100
    """
    bash_f.write("hmmer_wrapper.py -d outputCCS -i reads_of_insert.fasta --change-seqid --must-see-polyA --min-seqlen 100 -o reads_of_insert.53Aseen_trimmed_changeid.fa --cpus {cpus} || exit $?\n".format(cpus=args.cpus))
    bash_f.write("summarize_primer_info.py reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt > reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt.summary\n")

@check_file_or_run('reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa')    
def get_chimera(bash_f, args):
    bash_f.write("chimera_finder.py -d outputCCS_chimera --cpus {cpus} -i reads_of_insert.53Aseen_trimmed_changeid.fa || exit $?\n".format(cpus=args.cpus))
        
@check_file_or_run('reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa.blasr.sam')        
def run_BLASR(bash_f, args):
    bash_f.write("pbalign.py --maxHits 10 --hitPolicy all --algorithmOptions \"-nproc {cpus} -bestn 10 -nCandidates 10\"  reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa  {ref} reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa.blasr.sam\n".format(cpus=args.cpus, ref=args.ref))    

@check_file_or_run('reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff.log')
def run_gmap(bash_f, args):
    bash_f.write("gmap -D {gmap_db_dir} -d {gmap_db} -t {cpus} -f gff3_gene -n 0 reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa > reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff 2> reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff.log\n".format(cpus=args.cpus, gmap_db=args.gmap_db, gmap_db_dir=args.gmap_db_dir))
    

def main(args):
    f = open(args.cmd_filename, 'w')
    f.write("#!/bin/bash\n")
    f.write(". {0}/bin/activate\n".format(os.environ['VIRTUAL_ENV']))
    
    get_primer_info(f, args)
    get_chimera(f, args)
    run_BLASR(f, args)
    run_gmap(f, args)
    
    f.close()       
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--ref", default='/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/Gencode15/', help="Transcript reference location (default: Gencode15)")
    parser.add_argument("--gmap_db", default='hg19', help="GMAP db name (default: hg19)")
    parser.add_argument("--gmap_db_dir", default='/home/UNIXHOME/etseng/share/gmap_db/', help="GMAP db directory")
    parser.add_argument("--cpus", type=int, default=24, help="Number of CPUs (default: 24)")
    parser.add_argument("--force", action="store_true", default=False, help="Run all commands even if output files are present for some steps")
    parser.add_argument("--cmd_filename", default='cDNApipe.sh', help="Output CMD filename")
    
    args = parser.parse_args()
    main(args)     
