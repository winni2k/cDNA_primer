#!/usr/bin/env python
###########################################################################
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
###########################################################################
"""
Given an input_fasta file of non-full-length (partial) reads and
(unpolished) consensus isoforms sequences in ref_fasta, align reads to
consensus isoforms using BLASR, and then build up a mapping between
consensus isoforms and reads (i.e., assign reads to isoforms).
Finally, save
    {isoform_id: [read_ids],
     nohit: set(no_hit_read_ids)}
to an output pickle file.
"""


import os
import os.path as op
import sys
import logging
from cPickle import dump
from pbcore.io import FastaReader
from pbcore.util.Process import backticks
from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.Utils import realpath, touch
#from pbtools.pbtranscript.ice.ProbModel import ProbFromModel, ProbFromQV
from pbtools.pbtranscript.c_Prob import ProbFromModel, ProbFromQV
from pbtools.pbtranscript.ice.IceUtils import blasr_against_ref

#log = logging.getLogger(__name__)

def build_uc_from_partial(input_fasta, ref_fasta, out_pickle,
        sa_file=None, ccs_fofn=None, done_filename=None, blasr_nproc=12):
    """Align consensus isoforms in ref_fasta and reads in input_fasta,
    and save mappings between isoforms and reads to out_pickle.
    ccs_fofn --- If None, assume no quality value is available,
    otherwise, use QV from ccs_fofn.
    blasr_nproc --- equivalent to blasr -nproc, number of CPUs to use
    """
    input_fasta = realpath(input_fasta)
    m5_file = input_fasta + ".blasr"
    out_pickle = realpath(out_pickle)
    if sa_file is None:
        if op.exists(input_fasta + ".sa"):
            sa_file = input_fasta + ".sa"

    cmd = "blasr {i} ".format(i=input_fasta) + \
          "{r} -bestn 5 ".format(r=ref_fasta) + \
          "-nproc {n} -m 5 ".format(n=blasr_nproc) + \
          "-maxScore -1000 -minPctIdentity 85 -out {o} ".format(o=m5_file)
    if sa_file is not None and op.exists(sa_file):
        cmd += "-sa {sa}".format(sa=sa_file)

    logging.info("CMD: {cmd}".format(cmd=cmd))
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        errMsg = "Command failed: {cmd}\n{e}".format(cmd=cmd, e=_msg)
        logging.error(errMsg)
        raise RuntimeError(errMsg)

    if ccs_fofn is None:
        logging.info("Loading probability from model")
        probqv = ProbFromModel(.01, .07, .06)
    else:
        logging.info("Loading probability from QV in {f}".format(f=ccs_fofn))
        probqv = ProbFromQV(input_fofn=ccs_fofn, fasta_filename=input_fasta)

    logging.info("Calling blasr_against_ref ...")
    hitItems = blasr_against_ref(output_filename=m5_file,
                              is_FL=False,
                              sID_starts_with_c=True,
                              qver_get_func=probqv.get_smoothed,
                              ece_penalty=1,
                              ece_min_len=10,
                              same_strand_only=False)

    partial_uc = {}  # Maps each isoform (cluster) id to a list of reads
                     # which can map to the isoform
    seen = set()  # reads seen
    logging.info("Building uc from BLASR hits.")
    for h in hitItems:
        if h.ece_arr is not None:
            if h.cID not in partial_uc:
                partial_uc[h.cID] = []
            partial_uc[h.cID].append(h.qID)
            seen.add(h.qID)

    allhits = set(r.name.split()[0] for r in FastaReader(input_fasta))

    logging.info("Counting reads with no hit.")
    nohit = allhits.difference(seen)

    logging.info("Dumping uc to a pickle: {f}.".format(f=out_pickle))
    with open(out_pickle, 'w') as f:
        dump({'partial_uc': partial_uc, 'nohit': nohit}, f)

    os.remove(m5_file)

    done_filename = realpath(done_filename) if done_filename is not None \
            else out_pickle + '.DONE'
    logging.debug("Creating {f}.".format(f=done_filename))
    touch(done_filename)


def set_parser(parser):
    """Get arguments."""
    parser.add_argument("input_fasta",
                        help="Input fasta split file")
    parser.add_argument("ref_fasta",
                        help="Reference fasta file, most likely " +
                             "ref_consensus.fa from ICE output")
    parser.add_argument("out_pickle",
                        type=str,
                        default=None,
                        help="Output pickle file.")

    from pbtools.pbtranscript.PBTranscriptOptions import add_fofn_arguments
    parser = add_fofn_arguments(parser, ccs_fofn=True)

    parser.add_argument("--sa",
                        dest="sa_file",
                        default=None,
                        help="saffix array of ref_fasta")
    parser.add_argument("--done",
                        dest="done_filename",
                        help="An empty file generated to indicate that" +
                        "out_pickle is done.")
    parser.add_argument("--blasr_nproc",
                        dest="blasr_nproc",
                        type=int,
                        default=12,
                        help="blasr -nproc, number of CPUs [default: 12]")
    return parser


from pbcore.util.ToolRunner import PBToolRunner

class IcePartialRunner(PBToolRunner):
    """IcePartial Runner"""
    def __init__(self):
        desc = "Assign non-full-length reads in a splitted fasta file " + \
               "to unpolished isoform clusters"
        PBToolRunner.__init__(self, desc)
        set_parser(self.parser)

    def getVersion(self):
        """Get version string."""
        return get_version()

    def run(self):
        """Run"""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))

        args = self.args
        logging.info("Building uc from non-full-length reads.")
        build_uc_from_partial(input_fasta=args.input_fasta,
                              ref_fasta=args.ref_fasta,
                              out_pickle=args.out_pickle,
                              sa_file=args.sa_file,
                              ccs_fofn=args.ccs_fofn,
                              blasr_nproc=args.blasr_nproc)
        return 0


def main():
    """Main function"""
    runner = IcePartialRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
