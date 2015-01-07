import unittest
import os.path as op
from pbtools.pbtranscript.ice import IceUtils

class Test_ICEUtils(unittest.TestCase):
    """Test IceUtils."""
    def setUp(self):
        """Initialize."""
        self.outDir = op.join(op.dirname(op.dirname(op.abspath(__file__))),
                              "out")

    def test_sanity_check_gcon(self):
        """sanity_check_gcon."""
        self.assertTrue(IceUtils.gcon_py == "ice_gcontools.py")
        self.assertTrue(IceUtils.sanity_check_gcon() == IceUtils.gcon_py)

    #def test_sanity_check_sge(self):
    #    """sanity_check_sge."""
    #    self.assertTrue(IceUtils.sanity_check_sge(self.outDir))
