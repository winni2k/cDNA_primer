"""Test FastaSplitter"""

import unittest
import os.path as op
from pbcore.io import FastaReader
from pbtools.pbtranscript.io.FastaRandomReader import FastaRandomReader, \
        MetaSubreadFastaReader
from pbtools.pbtranscript.Utils import write_files_to_fofn
import hashlib

class TestFastaRandomReader(unittest.TestCase):
    """Class for testing FastaRandomReader."""
    def setUp(self):
        """Set up testDir, dataDir, outDir, stdoutDir"""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = self.rootDir
        self.dataDir = op.join(self.testDir, "data")
        self.outDir = op.join(self.testDir, "out")
        self.stdoutDir = op.join(self.testDir, "stdout")
        self.inFa = op.join(self.dataDir, "reads_of_insert.fasta")

    def testAll(self):
        """Test FastaRandomReader.keys() and __getitem__."""
        reads = [r for r in FastaReader(self.inFa)]
        names = [r.name for r in reads]
        seqs = [r.sequence for r in reads]

        frr = FastaRandomReader(self.inFa)
        self.assertTrue(set(frr.keys()) == set(names))

        self.assertTrue(False not in
                [frr[r.name].name == r.name for r in reads])

        self.assertTrue(False not in
                [frr[r.name].sequence == r.sequence for r in reads])

class TestMetaSubreadFastaReader(unittest.TestCase):
    """Class for testing MetaSubreadFastaReader."""
    def setUp(self):
        """Set up testDir, dataDir, outDir, stdoutDir"""
        self.testDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.dataDir = op.join(self.testDir, "data")
        self.outDir = op.join(self.testDir, "out")
        self.stdoutDir = op.join(self.testDir, "stdout")
        self.fa1 = op.join(self.dataDir, "test_meta_subreads_fasta_reader1.fasta")
        self.fa2 = op.join(self.dataDir, "test_meta_subreads_fasta_reader2.fasta")
        self.fofn = op.join(self.outDir, "test_meta_subreads_fasta_reader.fofn")

    def testAll(self):
        """Test FastaRandomReader.keys() and __getitem__."""
        write_files_to_fofn([self.fa1, self.fa2], self.fofn)
        reader = MetaSubreadFastaReader([self.fa1, self.fa2])
        subread_1 = "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/59/0_5071"
        subread_2 = "m130812_random_random_s1_p0/440/13280_16126"
        zmw_3 = "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/70"
        zmw_4 = "m130812_random_random_s1_p0/249"
        r1 = reader[subread_1][0]
        self.assertEqual(r1.name, subread_1)
        self.assertEqual(hashlib.md5(r1.sequence).hexdigest(), "8128261dd851ae285d029618739559e9")

        r2 = reader[subread_2][0]
        self.assertEqual(r2.name, subread_2)
        self.assertEqual(hashlib.md5(r2.sequence).hexdigest(), "451e5798a7f21cce80da27a03a8cb2c7")

        r3, r4 = reader[zmw_3]
        self.assertEqual(r3.name, "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/70/0_5538")
        self.assertEqual(r4.name, "m130812_185809_42141_c100533960310000001823079711101380_s1_p0/70/5587_5982")
        self.assertEqual(hashlib.md5(r3.sequence).hexdigest(), "4db2f6e35c83dd279a8f71a51ac50445")
        self.assertEqual(hashlib.md5(r4.sequence).hexdigest(), "1c1d080e9362a73ea2074f9a62fbd45e")

        r5 = reader[zmw_4][0]
        self.assertEqual(r5.name, "m130812_random_random_s1_p0/249/0_1339")
        self.assertEqual(hashlib.md5(r5.sequence).hexdigest(), "b20d3723a136aedc2f96f6f498ad3da0")



