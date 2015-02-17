Test pbtranscript.py
  $ CURDIR=$TESTDIR
  $ DATDIR=$CURDIR/../data
  $ OUTDIR=$CURDIR/../out
  $ STDDIR=$CURDIR/../stdout

#Test pbtranscript.py classify
  $ OUTFA=$OUTDIR/classify_out_1.fa 
  $ OUTCSV=$OUTDIR/classify_out_1.primer_info.csv 
  $ OUTSUMMARY=$OUTDIR/classify_out_1.classify_summary.txt 

  $ rm -f $OUTFA $OUTCSV $OUTSUMMARY
  $ pbtranscript.py classify $DATDIR/test.fa $OUTFA 
  $ echo $?
  0

  $ ls $OUTFA
  *classify_out_1.fa (glob)

  $ ls $OUTCSV
  *classify_out_1.primer_info.csv (glob)


  $ OUTFA=$OUTDIR/classify_out_2.fa 
  $ OUTCSV=$OUTDIR/classify_out_2.csv 
  $ OUTSUMMARY=$OUTDIR/classify_out_2.txt 
  $ rm -f $OUTFA $OUTCSV $OUTSUMMARY
  $ pbtranscript.py classify $DATDIR/reads_of_insert.fasta $OUTFA --report $OUTCSV --summary $OUTSUMMARY
  $ echo $?
  0

  $ ls $OUTFA
  *classify_out_2.fa (glob)

  $ ls $OUTCSV
  *classify_out_2.csv (glob)

  $ cat $OUTSUMMARY
  Number of reads of insert=22
  Number of five prime reads=11
  Number of three prime reads=13
  Number of poly-A reads=13
  Number of filtered short reads=0
  Number of non-full-length reads=13
  Number of full-length reads=9
  Number of full-length non-chimeric reads=9
  Average full-length non-chimeric read length=3951

  $ rm -f $OUTTXT $OUTFA
  $ INFA=$DATDIR/test_subset.fa 
  $ OUTFA=$OUTDIR/test_subset.fa 
  $ OUTTXT=$OUTDIR/test_subset.txt 
  $ pbtranscript.py subset $INFA $OUTFA --FL --nonChimeric
  $ echo $?
  0

  $ rm -f $OUTTXT
  $ pbtranscript.py subset $INFA $OUTTXT --FL --nonChimeric --printReadLengthOnly
  $ echo $?
  0
  $ cat $OUTTXT
  3889
  4194
  3950
  4409
  4106
  3784

  $ rm -f $OUTTXT
  $ pbtranscript.py subset $INFA $OUTTXT --chimeric --printReadLengthOnly
  $ echo $?
  0
  $ wc -l $OUTTXT | cut -f 1 -d ' ' 
  10
