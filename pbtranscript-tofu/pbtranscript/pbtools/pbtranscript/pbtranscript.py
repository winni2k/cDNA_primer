#!/usr/bin/env python
###############################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""This script defines class PBTranscript."""
import sys
import logging
from pbcore.util.ToolRunner import PBMultiToolRunner
from pbtools.pbtranscript.PBTranscriptException import PBTranscriptException
from pbtools.pbtranscript.Classifier import Classifier, ChimeraDetectionOptions
from pbtools.pbtranscript.ClusterOptions import IceOptions, SgeOptions, \
    IceQuiverHQLQOptions
from pbtools.pbtranscript.Cluster import Cluster
from pbtools.pbtranscript.SubsetExtractor import ReadsSubsetExtractor, \
    SubsetRules
from pbtools.pbtranscript.PBTranscriptOptions import \
    add_classify_arguments, add_cluster_arguments, add_subset_arguments
from pbtools.pbtranscript.__init__ import get_version


class PBTranscript(PBMultiToolRunner):

    """
    Class PBTranscript defines tool kits for cDNA analysis, including
    `classify` and `cluster`.
    `classify` - pbtranscript classifies reads from a fasta/q file.
    For each read, identify whether it is full length, whether 5', 3' and
    poly A tail have been found. The input fasta/q file is usually generated
    from RS_ReadsOfInsert protocol (e.g., reads_of_insert.fasta/q).
    `cluster` -  pbtranscript calls the ICE algorithm, which stands for
    'Iterarively Clustering and Error correction' to identify de novo
    consensus isoforms.
    """

    def __init__(self):
        desc = "Toolkit for cDNA analysis."
        super(PBTranscript, self).__init__(desc)
        subparsers = self.subParsers

        parser = subparsers.add_parser(
            'classify',
            description="Classify reads based on whether they are " +
                        "non-chimeric, full length and have their 5', " +
                        "3' and poly A tail seen.")
        # Add arguments for subcommand classify
        add_classify_arguments(parser)

        parser = subparsers.add_parser(
            'cluster',
            description='Discover consensus isoforms based on ' +
                        'quality controlled non-chimeric, ' +
                        'full length reads to reference genome.')
        # Add arguments for subcommand cluster
        add_cluster_arguments(parser, show_sge_env_name=True, show_sge_queue=True)

        parser = subparsers.add_parser(
            'subset',
            description='Subset annotated reads in Fasta format.')
        add_subset_arguments(parser)

    def getVersion(self):
        return get_version()

    def run(self):
        """Run classify, cluster, polish or subset."""
        cmd = self.args.subCommand
        try:
            if cmd == 'classify':
                opts = ChimeraDetectionOptions(
                    min_seq_len=self.args.min_seq_len,
                    min_score=self.args.min_score,
                    min_dist_from_end=self.args.min_dist_from_end,
                    max_adjacent_hit_dist=self.args.max_adjacent_hit_dist,
                    primer_search_window=self.args.primer_search_window,
                    detect_chimera_nfl=self.args.detect_chimera_nfl)

                obj = Classifier(reads_fn=self.args.readsFN,
                                 out_dir=self.args.outDir,
                                 out_reads_fn=self.args.outReadsFN,
                                 primer_fn=self.args.primerFN,
                                 primer_report_fn=self.args.primerReportFN,
                                 summary_fn=self.args.summary_fn,
                                 cpus=self.args.cpus,
                                 change_read_id=True,
                                 opts=opts,
                                 out_flnc_fn=self.args.flnc_fa,
                                 out_nfl_fn=self.args.nfl_fa,
                                 ignore_polyA=self.args.ignore_polyA,
                                 reuse_dom=self.args.reuse_dom)
                obj.run()
            elif cmd == 'cluster':
                ice_opts = IceOptions(quiver=self.args.quiver,
                                      use_finer_qv=self.args.use_finer_qv,
                                      targeted_isoseq=self.args.targeted_isoseq,
                                      ece_penalty=self.args.ece_penalty,
                                      ece_min_len=self.args.ece_min_len)
                sge_opts = SgeOptions(unique_id=self.args.unique_id,
                                      use_sge=self.args.use_sge,
                                      max_sge_jobs=self.args.max_sge_jobs,
                                      blasr_nproc=self.args.blasr_nproc,
                                      quiver_nproc=self.args.quiver_nproc,
                                      gcon_nproc=self.args.gcon_nproc,
                                      sge_env_name=self.args.sge_env_name,
                                      sge_queue=self.args.sge_queue)
                ipq_opts = IceQuiverHQLQOptions(qv_trim_5=self.args.qv_trim_5,
                                                qv_trim_3=self.args.qv_trim_3,
                                                hq_quiver_min_accuracy=self.args.hq_quiver_min_accuracy,
                                                hq_isoforms_fa=self.args.hq_isoforms_fa,
                                                hq_isoforms_fq=self.args.hq_isoforms_fq,
                                                lq_isoforms_fa=self.args.lq_isoforms_fa,
                                                lq_isoforms_fq=self.args.lq_isoforms_fq)

                obj = Cluster(root_dir=self.args.root_dir,
                              flnc_fa=self.args.flnc_fa,
                              nfl_fa=self.args.nfl_fa,
                              bas_fofn=self.args.bas_fofn,
                              ccs_fofn=self.args.ccs_fofn,
                              fasta_fofn=self.args.fasta_fofn,
                              out_fa=self.args.consensusFa,
                              sge_opts=sge_opts,
                              ice_opts=ice_opts,
                              ipq_opts=ipq_opts,
                              report_fn=self.args.report_fn,
                              summary_fn=self.args.summary_fn,
                              nfl_reads_per_split=self.args.nfl_reads_per_split)
                obj.run()

            elif cmd == 'subset':
                rules = SubsetRules(FL=self.args.FL,
                                    nonChimeric=self.args.nonChimeric)

                obj = ReadsSubsetExtractor(inFN=self.args.readsFN,
                                           outFN=self.args.outFN,
                                           rules=rules,
                                           ignore_polyA=self.args.ignore_polyA,
                                           printReadLengthOnly=self.args.printReadLengthOnly)
                obj.run()
            else:
                raise PBTranscriptException(cmd,
                                            "Unknown command passed to pbtranscript.py:"
                                            + self.args.subName)
        except Exception:
            logging.exception("Exiting pbtranscript with return code 1.")
            return 1
        return 0

if __name__ == "__main__":
    sys.exit(PBTranscript().start())
