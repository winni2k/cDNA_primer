"""Test pbtools.pbtranscript.Classifier."""

import unittest
import os
import os.path as op
from pbtools.pbtranscript.Classifier import Classifier, PBRead
from pbtools.pbtranscript.io.DOMIO import DOMRecord
from collections import namedtuple

import filecmp

class Test_Classifier(unittest.TestCase):
    """Test Classifier."""
    def setUp(self):
        """Set up test data."""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = op.join(self.rootDir, "")

    def test_processPrimers(self):
        """Test function _processPrimers()."""
        inPFN = op.join(self.testDir, "data/test_primers_in.fa")
        obj = Classifier()

        # Test on an artificial example.
        outPFN = op.join(self.testDir, "out/test_primers_out.fa")
        stdoutPFN = op.join(self.testDir, "stdout/test_primers_out.fa")
        obj._processPrimers(primer_fn=inPFN, window_size=50,
                            primer_out_fn=outPFN,
                            revcmp_primers=False)

        self.assertTrue(filecmp.cmp(outPFN, stdoutPFN))

        # Test on real PacBio primers.fa
        pbPFN = op.join(self.testDir, "data/primers.fa")

        # outPFN2 = primers.fa for primer detection.
        outPFN2 = op.join(self.testDir, "out/test_primers_out_2.fa")
        stdoutPFN2 = op.join(self.testDir, "stdout/test_primers_out_2.fa")
        obj._processPrimers(primer_fn=pbPFN, window_size=50,
                            primer_out_fn=outPFN2,
                            revcmp_primers=False)
        self.assertTrue(filecmp.cmp(outPFN2, stdoutPFN2))

        # outPFN3 = primers.fa for chimera detection.
        outPFN2 = op.join(self.testDir, "out/test_primers_out_3.fa")
        stdoutPFN2 = op.join(self.testDir, "stdout/test_primers_out_3.fa")
        obj._processPrimers(primer_fn=pbPFN, window_size=50,
                            primer_out_fn=outPFN2,
                            revcmp_primers=True)
        self.assertTrue(filecmp.cmp(outPFN2, stdoutPFN2))

    def test_chunkReads(self):
        """Test function _chunkReads(readsFN, chunkSize, chunkedReadsFNs)."""
        obj = Classifier()
        readsFN = op.join(self.testDir, "data/test_chunkReads_1.fa")
        chunkedReadsFN = op.join(self.testDir,
                                 "out/test_chunkReads_1.fa")
        if op.exists(chunkedReadsFN):
            os.remove(chunkedReadsFN)
        stdoutChunkedReadsFN = op.join(self.testDir,
                                       "stdout/test_chunkReads_1.fa")

        obj._chunkReads(readsFN, 10, [chunkedReadsFN])
        self.assertTrue(filecmp.cmp(chunkedReadsFN, stdoutChunkedReadsFN))

    def test_getBestFrontBackRecord(self):
        """Test function _parseBestFrontBackRecord()."""
        obj = Classifier()
        domFN = op.join(self.testDir, "data/test_parseHmmDom.dom")
        front, back = obj._getBestFrontBackRecord(domFN)
        # In the following, verify the front and back are equivalent
        # to stdout/test_parseHmmDom_dFront/Back.txt
        def prettystr(d):
            """Return Pretty print string for front & back."""
            return "\n".join(
                [key + ":\n" + "\n".join(
                    [k + ":" + str(v) for k, v in val.iteritems()])
                 for key, val in d.iteritems()])
        frontFN = op.join(self.testDir, "out/test_parseHmmDom_dFront.txt")
        backFN = op.join(self.testDir, "out/test_parseHmmDom_dBack.txt")
        f = open(frontFN, 'w')
        f.write(prettystr(front))
        f.close()
        f = open(backFN, 'w')
        f.write(prettystr(back))
        f.close()
        stdoutFrontFN = op.join(self.testDir,
                                "stdout/test_parseHmmDom_dFront.txt")
        stdoutBackFN = op.join(self.testDir,
                               "stdout/test_parseHmmDom_dBack.txt")
        self.assertTrue(filecmp.cmp(frontFN, stdoutFrontFN))
        self.assertTrue(filecmp.cmp(backFN, stdoutBackFN))

    def test_findPolyA(self):
        """Test function _findPolyA(seq, minANum, p3Start)."""
        obj = Classifier()
        seq1 = ("GTGAAGTAGGTGTCCCGCACCAAGGCACGGAGCCAGAGAGGTGTGGGTGC" +
                "TAAAAGCCACCCGTTAGGACCCAGAGCAGCTGAAGCTGGATGCGAAAGGA" +
                "TACAGGCTTAGTAGCCATGGAGACCAAACTGGAACAAATGCCGACTGGAA" +
                "AGTGTATCTTATAACTTATTAAATAAAATGTTTGCTCCACGAAAAAAAAA" +
                "AAAAAAAAAAAAAAGTACTCTGCGTTGATACCACTGCTT")
        seq2 = ("TGGTTGGTCGGCGTTTAGCTTTGTGAGGCTCCCTGAACAGAAACACTGTT" +
                "GGAAGAAGAGTCCCCTGACATCACCCAGCGTCAAGTGGGAGTTAGCCTCT" +
                "GAAGTTCAGTGTATCACGTTAATGCTAATATGCTTTGTGGTGGCAGAATT" +
                "TATTTTGGCTTTTTGTCATTTAGCCAAATTAAAGGCAAACGCGTTTCTAA" +
                "AAAAAAAAAAAAAAAAAAAAGTAGCTCTGCGTTTGATACCACTGCTT")
        seq3 = ("TATTTTGGCTTTTTGTCATTTAGCCAAATTAAAGGCAAACGCGTTTCTAA")
        self.assertEqual(obj._findPolyA(seq1), 188)
        self.assertEqual(obj._findPolyA(seq2), 196)
        self.assertEqual(obj._findPolyA(seq3), -1)

    def test_pickBestPrimerCombo(self):
        """Test funciton _pickBestPrimerCombo()."""
        obj = Classifier()
        domFN = op.join(self.testDir, "data/test_parseHmmDom.dom")
        front, back = obj._getBestFrontBackRecord(domFN)

        # Now pick up the best primer combo
        movie = "m131018_081703_42161_c100585152550000001823088404281404_s1_p0"
        rids = [movie + "/" + str(zmw) + "/ccs" for zmw in [43, 45, 54]]
        res = obj._pickBestPrimerCombo(
            front[rids[0]], back[rids[0]], [0, 1], 10)
        self.assertTrue(res[2] is None)
        self.assertTrue(res[3] is None)

        res = obj._pickBestPrimerCombo(
            front[rids[1]], back[rids[1]], [0, 1], 10)

        fw = DOMRecord("F1", movie + "/45/ccs", 33.0, 0, 30, 31, 0, 30, 100)
        rc = DOMRecord("R1", movie + "/45/ccs", 27.2, 0, 25, 25, 0, 25, 100)
        self.assertEqual(res[0], 1)
        self.assertEqual(res[1], "+")
        self.assertTrue(str(fw) == str(res[2]))
        self.assertTrue(str(rc) == str(res[3]))

        res = obj._pickBestPrimerCombo(
            front[rids[2]], back[rids[2]], [0, 1], 10)
        rc = DOMRecord("R1", movie + "/54/ccs", 22.3, 0, 25, 25, 0, 27, 100)
        self.assertEqual(res[0], 1)
        self.assertEqual(res[1], "+")
        self.assertTrue(res[2] is None)
        self.assertTrue(str(res[3]) == str(rc))

    def test_PBRead(self):
        """Test class PBRead."""
        A = namedtuple('A', 'name sequence')
        x = A("movie/10/2_100", "A" * 98)
        y = PBRead(x)
        self.assertEqual((y.movie, y.zmw, y.start, y.end, y.isCCS),
                ("movie", 10, 2, 100, False))
        x = A("movie/10", '')
        y = PBRead(x)
        self.assertEqual((y.movie, y.zmw, y.isCCS),
                ("movie", 10, True))
        x = A("movie/10/ccs", '')
        y = PBRead(x)
        self.assertEqual((y.movie, y.zmw, y.isCCS),
                ("movie", 10, True))

if __name__ == "__main__":
    unittest.main()
