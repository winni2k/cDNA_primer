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

import os, sys
from collections import defaultdict
import numpy as np
import pbtools.pbtranscript.c_basQV as c_basQV
from pbcore.io import BasH5Reader

class h5_wrapper:
    """
    Wrap the .1.ccs.h5, .2.ccs.h5, .3.ccs.h5 to give the illusion of one
    (like what .bas.h5 does for .bax.h5)
    """
    def __init__(self, file_prefix, suffix='.ccs.h5'):
        """
        Expects <prefix>.1.ccs.h5, .2.ccs.h5, and .3.ccs.h5
        or .1.bax.h5, .2.bax.h5, .3.bax.h5
        """
        assert suffix in ('.ccs.h5', '.bax.h5')
        self.files = [file_prefix+'.1'+suffix, file_prefix+'.2'+suffix, file_prefix+'.3'+suffix]
        assert all(os.path.exists(file) for file in self.files)

        self.hn_range = [(0,54493),(54494,108987),(108988,170000)]
        #self.get_hn_range()

    def get_hn_range(self):
        for file in self.files:
            print >> sys.stderr, "getting holeNumber range for", file
            bas = BasH5Reader(file)
            _lo = bas.sequencingZmws[0]
            _hi = bas.sequencingZmws[-1]
            self.hn_range.append((_lo, _hi))

    def __getitem__(self, seqid):
        if seqid.count('/') == 1: hn = int(seqid.split('/')[1])
        elif seqid.count('/') == 2: hn = int(seqid.split('/')[1])
        else: raise Exception, "Cannot recognize", seqid
        for i in xrange(3):
            if self.hn_range[i][0] <= hn <= self.hn_range[i][1]: return self.files[i]
        raise Exception, "Unlocated holeNumber", hn


class basQVcacher:
    qv_names = ['InsertionQV', 'SubstitutionQV', 'DeletionQV']
    def __init__(self):
        self.bas_dict = {}
        self.bas_files = {}
        self.qv = {} # subread seqid --> qv_name --> list of qv (transformed to prob)
        self.window_size = None # smoothing window size, set when presmooth() is called

    def get(self, seqid, qv_type, position=None):
        if position is None:
            return self.qv[seqid][qv_type]
        else:
            return self.qv[seqid][qv_type][position]

    def get_smoothed(self, seqid, qv_type, position=None):
        if position is None:
            return self.qv[seqid][qv_type+'_smoothed']
        else:
            return self.qv[seqid][qv_type+'_smoothed'][position]

    def add_bash5(self, bash5_filename):
        basename = os.path.basename(bash5_filename)        
        if bash5_filename.endswith('.bax.h5'):
            movie = basename[:-9]
            if movie not in self.bas_files:
                self.bas_files[movie] = h5_wrapper(bash5_filename[:-9], suffix='.bax.h5')
        elif bash5_filename.endswith('.1.ccs.h5') or bash5_filename.endswith('.2.ccs.h5') or bash5_filename.endswith('.3.ccs.h5'):
            movie = basename[:-9]
            if movie not in self.bas_files:
                self.bas_files[movie] = h5_wrapper(bash5_filename[:-9])
        elif bash5_filename.endswith('.ccs.h5'): # a single .ccs.h5 (post 150k runs), treat the same as .bas.h5
            movie = basename[:-7]
            self.bas_files[movie] = defaultdict(lambda: bash5_filename)
        else:
            assert bash5_filename.endswith('.bas.h5') 
            movie = basename[:-7]
            self.bas_files[movie] = defaultdict(lambda: bash5_filename)

    def precache(self, seqids):
        """
        """
        # for subread ex: m120407_063017_42141_c100320212550000001523017409061204_s1_p0/13/2571_3282
        # for CCS ex:  m120407_063017_42141_c100320212550000001523017409061204_s1_p0/13/300_10_CCS
        
        # split seq IDs by bas filename
        bas_job_dict = defaultdict(lambda: [])  # bas file --> list of seqids belonging to that bas file
        for seqid in seqids:
            movie, hn, s_e = seqid.split('/')
            hn = int(hn)
            bas_file = self.bas_files[movie][seqid]
            bas_job_dict[bas_file].append(seqid)

        for bas_file, seqids in bas_job_dict.iteritems():
            c_basQV.precache_helper(bas_file, seqids, basQVcacher.qv_names, self.qv)

    def presmooth(self, seqids, window_size):
        """
        precache MUST BE already called! Otherwise will have error!
        """
        self.window_size = window_size
        for seqid in seqids:
            for qv_name in basQVcacher.qv_names:
                self.qv[seqid][qv_name+'_smoothed'] = c_basQV.maxval_per_window(self.qv[seqid][qv_name], window_size)

    def remove_unsmoothed(self):
        for k,v in self.qv.iteritems():
            for qv_name in basQVcacher.qv_names:
                try:
                    del self.qv[k][qv_name]
                except KeyError:
                    pass # may have already been deleted. OK.
                
                


