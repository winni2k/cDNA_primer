These scripts are for identifying potential full-length (FL) subreads or CCS reads using the 5' and 3' primer ligated to the transcripts during the cDNA library preparation step.

IMPORTANT: usage of the scripts is detailed in the wiki_ section. Please read it!!

.. _wiki: https://github.com/Magdoll/cDNA_primer/wiki


======================================================================                    
Identifying full-length subreads/CCS reads using cDNA kit primers
======================================================================

See this page_ on how to use the full-length identification scripts. 

.. _page: https://github.com/Magdoll/cDNA_primer/wiki/How-to-identify-full-length-transcripts-in-PacBio-data


usage: hmmer_wrapper.py
       [-h] [-p PRIMER_FILENAME] [-i INPUT_FILENAME] [-d DIRECTORY]
       [-k PRIMER_SEARCH_WINDOW] [--cpus CPUS] [--left-nosee-ok]
       [--right-nosee-ok] [--output-anyway] [--change-seqid]
       [--min-seqlen MIN_SEQLEN] [--min-score MIN_SCORE] -o OUTPUT_FILENAME

 This script requires phmmer from HMMER 3.0.
 If the output directory already exists, will skip running phmmer and directory go to primer trimming.
 If you want to re-run HMMER you must first delete the output directory manually.
 Refer to wiki: https://github.com/PacificBiosciences/cDNA_primer/wiki for more details.

optional arguments:
  -h, --help            show this help message and exit

HMMER options:
  -p PRIMER_FILENAME, --primer_filename PRIMER_FILENAME
                        Primer fasta file
  -i INPUT_FILENAME, --input_filename INPUT_FILENAME
                        Input fasta file (usually filtered_subreads.fasta or filtered_CCS_subreads.fasta)
  -d DIRECTORY, --directory DIRECTORY
                        Directory to store HMMER output (default: output/)
  -k PRIMER_SEARCH_WINDOW, --primer_search_window PRIMER_SEARCH_WINDOW
                        Search in the first/last k-bp for primers. Must be longer than the longest primer. (default: 100)
  --cpus CPUS           Number of CPUs to run HMMER (default: 8)

Primer trimming options:
  --left-nosee-ok       OK if 5' end not detected (default: off)
  --right-nosee-ok      OK if 3' end not detected (default: off)
  --output-anyway       Still output seqs w/ no primer (default: off)
  --change-seqid        Change seq id to reflect trimming (default: off)
  --min-seqlen MIN_SEQLEN
                        Minimum seqlength to output (default: 50)
  --min-score MIN_SCORE
                        Minimum bit score for primer hit (default: 10)
  -o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
                        Output fasta filename




===========================================================                    
Extra filtering to eliminate subreads with missed adapters
===========================================================

If SMRTbell adapters are missed, sometimes it'll still be considered full-length by barcode_trimmer.py (especially
when the 5' and 3' primers are identical or highly similar). To further eliminate these subreads, after running
*hmmer_wrapper.py*, run this on the remaining FL reads::
    chimera_finder.py -d <output_dir> --cpus <cpus> -i <FL fasta filename>




[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/7a0be18cf32fa6aca5c388e1a18cdf6d "githalytics.com")](http://githalytics.com/PacificBiosciences/cDNA_primer)


