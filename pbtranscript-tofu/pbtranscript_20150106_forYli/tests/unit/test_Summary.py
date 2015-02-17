"""Test pbtools.pbtranscript.Classifier."""
import unittest
import os.path as op
from pbtools.pbtranscript.io.Summary import ClassifySummary, ClusterSummary

import filecmp

class Test_ClassifySummary(unittest.TestCase):
    """Test ClassifySummary."""
    def setUp(self):
        """Set up test data."""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = op.join(self.rootDir, "")

    def test_write(self):
        """Test ClassifySummary.write."""
        outFN = op.join(self.testDir, "out/test_ClassifySummary.txt")
        stdoutFN = op.join(self.testDir, "stdout/test_ClassifySummary.txt")

        obj = ClassifySummary()
        obj.num_reads = 100
        obj.num_5_seen = 90
        obj.num_3_seen = 70
        obj.num_polyA_seen = 70
        obj.num_filtered_short_reads = 10
        obj.num_nfl = 50
        obj.num_fl = 40
        obj.num_flnc = 39
        obj.num_flc = 1
        obj.num_flnc_bases = 10001

        obj.write(outFN)
        self.assertTrue(filecmp.cmp(outFN, stdoutFN))

class Test_ClusterSummary(unittest.TestCase):
    """Test ClusterSummary."""
    def setUp(self):
        """Set up test data"""
        self.testDir = op.dirname(op.dirname(op.abspath(__file__)))

    def test_write(self):
        """Test ClusterSummary.write."""
        outFN = op.join(self.testDir, "out/test_ClusterSummary.txt")
        stdoutFN = op.join(self.testDir, "stdout/test_ClusterSummary.txt")

        obj = ClusterSummary()
        obj.numConsensusIsoforms = 97
        obj.numTotalBases = 97 * 3945

        obj.write(outFN)
        self.assertTrue(filecmp.cmp(outFN, stdoutFN))

