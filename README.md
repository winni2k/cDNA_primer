Last Updated: 01/31/2014

## About This Repository

The scripts in this repository is developed for the purpose of analyzing transcriptome data generated using the PacBio [Iso-Seq(TM) protocol](http://www.smrtcommunity.com/Share/Protocol?id=a1q70000000HqSvAAK&strRecordTypeName=Protocol). 


This code is not part of the official PacBio software package and is developed solely by Elizabeth Tseng. Plans to incorporate this code into the official software pipeline is underway and is expected to be available in 2.2 release.


The code is divided into two sections:

*Quality Control* --- scripts for (a) identifying full-length reads and (b) identifying artificial chimeric reads. The scripts are in the [scripts/](https://github.com/PacificBiosciences/cDNA_primer/tree/master/scripts) directory and tutorial #1-#5 from the [wiki](https://github.com/PacificBiosciences/cDNA_primer/wiki) covers its usage. Read scripts/INSTALL for how to compile code.

*ICE* --- scripts for (a) isoform-level clustering for consensus calling; (b) recruitment of non-FL reads and calling [Quiver](https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/HowToQuiver.rst) for final consensus polishing. ICE is released in a [separate tarball](https://github.com/PacificBiosciences/cDNA_primer/releases) and ICE tutorial #1-#2 from the [wiki](https://github.com/PacificBiosciences/cDNA_primer/wiki) covers its usage. The QC scripts are a pre-requisite for running ICE.


## About the Iso-Seq(TM) Application

Iso-Seq (Isoform Sequencing) Applicatino is the PacBio trademark name referring to its developed applications for transcriptome sequencing. Please refer to the resources below as well as the wiki for more information:

* [Iso-Seq Webinar](https://s3.amazonaws.com/files.pacb.com/Customer+Webinars/MCF-7+Transcriptome+Iso-Seq+Webinar+01+22+14.wmv) (recorded 1/22/2014)
* [Iso-Seq Webinar Slides](https://s3.amazonaws.com/files.pacb.com/pdf/Iso-Seq+Bioinformatics+Analysis+of+the+Human+MCF-7+Transcriptome.pdf)
* [Iso-Seq Webinar Q&A](https://s3.amazonaws.com/files.pacb.com/Customer+Webinars/Iso-Seq+Webinar+Q%26A.pdf)


## License

Standard PacBio Open Source License that is associated with this package:

```

#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$
```
