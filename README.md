[DevNet](https://github.com/PacificBiosciences/cDNA_primer/wiki) | <a href="mailto:devnet@pacificbiosciences.com">Contact Us</a> | [Terms of Use](http://pacbiodevnet.com/Terms_of_Use.html) | [Trademarks](http://pacb.com/terms-of-use/index.html#trademarks)

Last Updated: 05/12/2014



## Latest News

The latest *pbtranscript-tofu* version: 0.2.tofu.134691

The scripts in this repository are now compatible with SMRTAnalysis 2.2. It is a "beta" version of the RS_IsoSeq protocol that is supported in 2.2, meaning that it extends the existing RS_IsoSeq code base, and is an unofficial, developemental version of RS_IsoSeq. This beta version is only accessible via command line. The command line name of the official RS_IsoSeq is *pbtranscript*, to differentiate from the official version, in the tutorials we will refer to the beta version provided on this site as *pbtranscript-tofu*. 

Functionalities provided in *pbtranscript-tofu* are highly developmental and may not make it into the next version of official RS_IsoSeq.

Latest version of *pbtranscript-tofu* contains the following extensions:

* Support for identifying full-length transcripts that does not have polyA tail (ex: RT-PCR transcripts)
* Faster I/O for base QV reading in clustering, resulting in overall speedup
* Minor changes to better support SMRTAnalysis 2.2 framework


## About This Repository

The scripts in this repository is developed for the purpose of analyzing transcriptome data generated using the PacBio(R) [Iso-Seq(TM) protocol](http://www.smrtcommunity.com/Share/Protocol?id=a1q70000000HqSvAAK&strRecordTypeName=Protocol). 


This code is not part of the official PacBio software package and is developed solely by Elizabeth Tseng. The official software version is the RS_IsoSeq protocol. Use this code at your own risk.


## About Iso-Seq(TM)

The Iso-Seq (Isoform Sequencing) protocol refers to PacBio’s proprietary methods and applications for transcriptome sequencing. Please refer to the resources below, as well as the wiki for more information:

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
