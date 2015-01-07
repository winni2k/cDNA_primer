"""Test pbtools.pbtranscript.SubsetExtractor."""

import unittest
import os.path as op
from pbtools.pbtranscript.SubsetExtractor import SubsetRules, \
    ReadsSubsetExtractor
from pbtools.pbtranscript.io.ReadAnnotation import ReadAnnotation
from pbcore.io import FastaReader
import filecmp

class Test_SubsetExtractor(unittest.TestCase):
    """Test SubsetExtractor."""
    def setUp(self):
        """Set up test data."""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = op.join(self.rootDir, "")

    def test_satisfy(self):
        """Test function satisfy()."""
        inFN = op.join(self.testDir, "data/test_subset.fa")
        reads = []
        with FastaReader(inFN) as reader:
            reads = [x for x in reader]

        rules = SubsetRules(1, 1) # Full-length, non-chimeric
        obj = ReadsSubsetExtractor("in", "out", rules, True)

        ans = [ReadAnnotation.fromString(r.name) for r in reads]
        res = [obj.satisfy(an, rules) for an in ans]
        expected = [1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1]
        self.assertTrue(res == expected)

    def test_run(self):
        """Test function run()."""
        inFN = op.join(self.testDir, "data/test_subset.fa")
        outFN = op.join(self.testDir, "out/test_subset_unit.fa")
        stdoutFN = op.join(self.testDir, "stdout/test_subset_unit.fa")

        rules = SubsetRules(1, 1) # Full-length, non-chimeric
        obj = ReadsSubsetExtractor(inFN, outFN, rules, True)
        obj.run()
        self.assertTrue(filecmp.cmp(outFN, stdoutFN))





