#!/usr/bin/env python

import unittest
import os.path as op
import filecmp
from pbtools.pbtranscript.Utils import cat_files, filter_sam

class TestUtils(unittest.TestCase):
    """Test pbtools.pbtranscript.Utils"""
    def setUp(self):
        """Initialize."""
        self.root_dir = op.dirname(op.dirname(op.abspath(__file__)))
        self.data_dir = op.join(self.root_dir, "data")
        self.out_dir = op.join(self.root_dir, "out")
        self.stdout_dir = op.join(self.root_dir, "stdout")

    def test_cat_files(self):
        """Test cat_files."""
        fn_1 = op.join(self.data_dir, "primers.fa")
        fn_2 = op.join(self.data_dir, "test_phmmer.fa")
        out_fn_1 = op.join(self.out_dir, "test_cat_1")
        out_fn_2 = op.join(self.out_dir, "test_cat_2")

        std_out_fn_2 = op.join(self.stdout_dir, "test_cat_2")

        cat_files(src=[fn_1], dst=out_fn_1)
        cat_files(src=[fn_1, fn_2], dst=out_fn_2)
        self.assertTrue(filecmp.cmp(out_fn_1, fn_1))
        self.assertTrue(filecmp.cmp(out_fn_2, std_out_fn_2))

    def test_filter_sam(self):
        """Test filter_sam."""
        in_sam = op.join(self.data_dir, "test_filter_sam.sam")
        out_sam = op.join(self.out_dir, "test_filter_sam.sam")
        stdout_sam = op.join(self.stdout_dir, "test_filter_sam.sam")
        filter_sam(in_sam=in_sam, out_sam=out_sam)
        self.assertTrue(filecmp.cmp(out_sam, stdout_sam))
