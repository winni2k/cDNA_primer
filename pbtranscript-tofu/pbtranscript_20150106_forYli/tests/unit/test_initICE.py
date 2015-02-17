"""Test initICE."""

import unittest

class TestInitICE(unittest.TestCase):
    """Class for testing initICE."""
    def setUp(self):
        """Set up testDir, dataDir, outDir, stdoutDir."""
        self.rootDir = op.dirname(op.dirname(op.abspath(__file__)))
        self.testDir = op.join(self.rootDir, "")
