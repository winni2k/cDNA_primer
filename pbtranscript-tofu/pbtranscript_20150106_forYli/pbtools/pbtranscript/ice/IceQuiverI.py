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

"""
Overview:
    pbtranscript cluster contains two main components:
    * (1) ICE (iterative clustering and error correction) to predict
      unpolished consensus isoforms.
    * (2) Polish, to use nfl reads and quiver to polish those predicted
      unpolished isoforms. Polish contains three steps:
      + (2.1) IceAllPartials (ice_partial.py all)
              Align and assign nfl reads to unpolished isoforms, and
              save results to a pickle file.
      + (2.2) IceQuiver (ice_quiver.py i and ice_quiver.py merge)
              Call quiver to polish each isoform based on alignments
              created by mapping its associated fl and nfl reads to
              this isoform.
      + (2.3) IceQuiverPostprocess (ice_quiver.py postprocess)
              Collect and post process quiver results, and classify
              HQ/LQ isoforms.

    In order to handle subtasks by SMRTPipe instead of pbtranscript
    itself, we will refactor the polish phase including
    (2.1) (2.2) and (2.3). The refactor of (2.1) is described in
    ice_partial.py.

    (2.2) IceQuiver will be refactored to
       + (2.2.1) IceQuiverI (ice_quiver.py i)
                 Split all unpolished isoforms into N chunks and
                 call Quiver to polish isoforms of the i-th chunk
                 at a time
       + (2.2.2) IceQuiverMerge (ice_quiver.py merge)
                 When all splitted quiver jobs are done,
                 collect all submitted jobs and save to
                 root_dir/log/submitted_quiver_jobs.txt
    (2.3) IceQuiverPostProcess's entry will be renamed from
          ice_post_quiver.py to:
       + (2.3.1) ice_quiver.py postprocess

    e.g., IceQuiver (2.2) =
               ice_quiver_i.py root_dir 0 N   --> process the 0 th chunk
            +  ...
            +  ice_quiver_i.py root_dir N-1 N --> process the N-1 th chunk
            +  ice_quiver_merge.py root_dir N --> merge all polished isoforms
                                                  from N chunks

   *** Here we are focusing on (2.2.1) 'ice_quiver.py i' ***

Description:
    (2.2.2) ice_quiver.py i

    Assumption:
     * Phase (1) ICE is done, unpolished isoforms are created, fl reads
       are assigned to isoforms, and saved to a pickle (i.e., final.pickle)
     * Step (2.1) IceAllPartials is done, all nfl reads are assigned
       to unpolished isoforms, saved to a pickle (i.e., nfl_all_pickle_fn)

    Process:
       Given root_dir, i, and N, process the (i / N)-th workload of
       quiver jobs.

       How to divide all quiver jobs into N workloads?
       For all clusters (unpolished consensus isoforms), sort clusters by
       ids and then divide them evenly into N parts.

    Input:
        Positional:
            root_dir, an output directory for running pbtranscript cluster.
            i, the i-th chunk of quiver workload
            N, the total number of quiver workload chunks.
        Nont-positional:
            bas_fofn, fofn of input bas/bax.h5 files, required.
            fasta_fofn, fofn of fasta files extracted from bas/bax.h5 files, required.

    Output:
        All unpolished consensus isoforms in the (i / N) chunk will be polished
        by Quiver.

    hierarchy:
        pbtranscript = iceiterative

        pbtranscript --quiver = iceiterative + \
                                ice_polish.py

        ice_polish.py =  ice_make_fasta_fofn.py + \
                         ice_partial.py all + \
                         ice_quiver.py all

        ice_partial.py all = ice_partial.py split + \
                             ice_partial.py i + \
                             ice_partial.py merge

        (ice_partial.py one --> only apply ice_partial on a given input fasta)

        ice_quiver.py all = ice_quiver.py i + \
                            ice_quiver.py merge + \
                            ice_quiver.py postprocess

    Example:
        ice_quiver.py i root_dir {i} N --bas_fofn=bas_fofn \
                        --fasta_fofn=fasta_fofn

"""

from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.PBTranscriptOptions import add_fofn_arguments, \
    add_sge_arguments, add_cluster_root_dir_as_positional_argument
from pbtools.pbtranscript.ice.IceQuiver import IceQuiver


def add_ice_quiver_i_arguments(parser):
    """Add IceQuiverI (ice_quiver.py i) arguments."""
    parser = add_cluster_root_dir_as_positional_argument(parser)

    helpstr = "Call quiver to polish clusters in the (i / N)-th chunk."
    parser.add_argument("i", help=helpstr, type=int)

    helpstr = "Divide clusters into N chunks."
    parser.add_argument("N", help=helpstr, type=int)

    parser = add_fofn_arguments(parser, bas_fofn=True, fasta_fofn=True)
    parser = add_sge_arguments(parser, quiver_nproc=True, blasr_nproc=True)
    return parser


class IceQuiverI(object):

    """Divide all unpolished consensus isoforms into N chunks,
    and call quiver to polish isoforms of the i-th chunk."""

    desc = "Divide all predicted consensus isoforms into N chunks, " + \
           "and call quiver to polish isoforms of the i-th chunk."

    prog = "ice_quiver.py i "

    def __init__(self, root_dir, i, N, bas_fofn, fasta_fofn, sge_opts):
        self.root_dir = root_dir
        self.i = int(i)
        self.N = int(N)
        self.bas_fofn = bas_fofn
        self.fasta_fofn = fasta_fofn
        self.sge_opts = sge_opts
        if bas_fofn is None or fasta_fofn is None:
            raise ValueError("Both bas_fofn and fasta_fofn are required.")

    def getVersion(self):
        """Return version string."""
        return get_version()

    def cmd_str(self):
        """Return a cmd string of IceQuiverI (ice_quiver.py i)."""
        return self._cmd_str(root_dir=self.root_dir, i=self.i, N=self.N,
                             bas_fofn=self.bas_fofn, fasta_fofn=self.fasta_fofn,
                             sge_opts=self.sge_opts)

    def _cmd_str(self, root_dir, i, N, bas_fofn, fasta_fofn, sge_opts):
        """Return a cmd string of IceQuiverI (ice_quiver.py i)."""
        cmd = self.prog + \
            "{d} ".format(d=root_dir) + \
            "{i} ".format(i=i) + \
            "{N} ".format(N=N) + \
            "--bas_fofn={bas_fofn} ".format(bas_fofn=bas_fofn) + \
            "--fasta_fofn={fasta_fofn} ".format(fasta_fofn=fasta_fofn)
        cmd += sge_opts.cmd_str(show_blasr_nproc=True, show_quiver_nproc=True)
        return cmd

    def run(self):
        """Run"""
        cmd_str = self.cmd_str()

        iceq = IceQuiver(root_dir=self.root_dir,
                         bas_fofn=self.bas_fofn,
                         fasta_fofn=self.fasta_fofn,
                         sge_opts=self.sge_opts,
                         prog_name="ice_quiver_{i}of{N}".
                                   format(i=self.i, N=self.N))

        iceq.add_log(cmd_str)
        iceq.add_log("root_dir: {d}.".format(d=self.root_dir))
        iceq.add_log("Index of chunks: i = {i}.".format(i=self.i))
        iceq.add_log("Total number of chunks: N = {N}.".format(N=self.N))

        iceq.validate_inputs()

        iceq.process_chunk_i(i=self.i, num_chunks=self.N)
        iceq.close_log()


# import logging
# import sys
# import os.path as op
# from pbcore.util.ToolRunner import PBToolRunner
# from pbtools.pbtranscript.ClusterOptions import SgeOptions
#
#
# class IceQuiverIRunner(PBToolRunner):
#
#     """IceQuiverI runner"""
#
#     def __init__(self):
#         PBToolRunner.__init__(self, IceQuiverI.desc)
#         add_ice_quiver_i_arguments(self.parser)
#
#     def getVersion(self):
#         """Return version string."""
#         return get_version()
#
#     def run(self):
#         """Run"""
#         logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
#                                                 v=get_version()))
#         cmd_str = ""
#         try:
#             args = self.args
#             sge_opts = SgeOptions(unique_id=args.unique_id,
#                                   use_sge=args.use_sge,
#                                   max_sge_jobs=args.max_sge_jobs,
#                                   blasr_nproc=args.blasr_nproc,
#                                   quiver_nproc=args.quiver_nproc)
#
#             iceqi = IceQuiverI(root_dir=args.root_dir,
#                                bas_fofn=args.bas_fofn,
#                                fasta_fofn=args.fasta_fofn,
#                                sge_opts=sge_opts,
#                                prog_name="ice_quiver_{i}of{N}".
#                                format(i=args.i, N=args.N))
#
#             cmd_str = iceqi.cmd_str()
#             iceqi.run()
#         except:
#             logging.exception("Exiting {cmd_str} with return code 1".
#                               format(cmd_str=cmd_str))
#             return 1
#         return 0
#
#
# def main():
#     """Main function."""
#     runner = IceQuiverIRunner()
#     return runner.start()
#
# if __name__ == "__main__":
#     sys.exit(main())
