#!/usr/bin/env python
import os, sys, subprocess

species = sys.argv[1]

# ----------------------- SETTINGS --------------------- #
NUM_CPUS = 12
gmap_db_dir = '/home/UNIXHOME/etseng/share/gmap_db'

# option for GMAP and BLASR
if species == 'hg19':
    transcript_ref_dir = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/Gencode15/'
    transcript_ref_fasta = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/Gencode15/sequence/Gencode15.fasta'
    gmap_db_name = 'hg19'
elif species == 'rn5':
    transcript_ref_dir = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/rat_UCSC'
    transcript_ref_fasta = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/rat_UCSC/sequence/rat_UCSC.fasta'
    gmap_db_name = 'rn5'
elif species == 'mm10':
    transcript_ref_dir = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/mouse_UCSC'
    transcript_ref_fasta = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/mouse_UCSC/sequence/mouse_UCSC.fasta'
    gmap_db_name = 'mm10'
else:
    print >> sys.stderr, "species not specified or unknown! quit!"
    sys.exit(-1)

cmd_cDNA = "generate_cDNApipe_bash.py --ref {0} --gmap_db {1} --gmap_db_dir {2} --cpus {3} --cmd_filename cDNApipe.sh".format(transcript_ref_dir, gmap_db_name, gmap_db_dir, NUM_CPUS)
cmd_eval = "eval_cDNApipe_results.py --ref_fasta_filename {0} > evaled_summary.txt".format(transcript_ref_fasta)

# ----------------------- SETTINGS --------------------- #

# SANITY CHECK for runs/, smrtpipe/, primer_match
if not os.path.exists('runs/'):
    print >> sys.stderr, "runs/ directory does not exist. Quit!"
    sys.exit(-1)
if not os.path.exists('smrtpipe/'):
    print >> sys.stderr, "smrtpipe/ directory does not exist. Quit!"
    sys.exit(-1)
if not os.path.exists('primer_match/'):
    print >> sys.stderr, "primer_match/ directory does not exist. Quit!"
    sys.exit(-1)
if not os.path.exists('primer_match/primers.fa'):
    print >> sys.stderr, "primer_match/primers.fa does not exist. Quit!"
    sys.exit(-1)

# FOR EACH <name> in runs/
#  find the corressponding smrtpipe/<name>
#  create primer_match/<name>
#  symbolically link reads_of_insert.fasta to primer_match/<name>
#  generate cDNA script
#  submit cDNA script
for name in os.listdir('runs/'):
    smrt_d = os.path.join('smrtpipe', name, 'data')
    fa_ccs1 = os.path.join(smrt_d, 'reads_of_insert.fasta')
    if not os.path.exists(fa_ccs1):
        print >> sys.stderr, "{0} does not yet exist. Skipping {1}".format(fa_ccs1, name)
        continue

    pm_d = os.path.join('primer_match', name)
    if not os.path.exists(pm_d): os.makedirs(pm_d)
    else: continue # skipping pm_d because already exists
    fa_primer = os.path.join(pm_d, 'primers.fa')
    if not os.path.exists(fa_primer):
        if os.path.exists('primer_match/primers.fa'):
            cmd = "cp primer_match/primers.fa " + fa_primer
            if os.system(cmd)!=0:
                print >> sys.stderr, "Trouble with {0}.".format(cmd)
                continue
        else:
            print >> sys.stderr, "Need primer file {0}!".format(fa_primer)
            continue

    fa_ccs2 = os.path.join(pm_d, 'reads_of_insert.fasta')
    if not os.path.exists(fa_ccs2) and os.system("ln -s " + os.path.abspath(fa_ccs1) + " " + fa_ccs2)!=0:
        print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(fa_ccs1, name)
        continue


    cwd = os.popen("pwd").read().strip()
    os.chdir(pm_d)
    subprocess.check_call(cmd_cDNA, shell=True)
    with open('cDNApipe.sh', 'a') as f: 
        f.write(cmd_eval + '\n')
    cmd = "qsub -cwd -S /bin/bash -pe smp {cpus} cDNApipe.sh".format(cpus=NUM_CPUS)
    print >> sys.stderr, "submitting job for ", pm_d
    subprocess.check_call(cmd, shell=True)
    os.chdir(cwd)
