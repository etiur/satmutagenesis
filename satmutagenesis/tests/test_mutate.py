"""
This module provides test of the mutate_pdb.py's functions
"""

class TestMutagenesis:
    """
    Test the class Mutagenesis on mutate.py
    """
    def test_mutate(self):
        """
        It checks if the pmx library is working correctly
        """
        from mutate_pdb import Mutagenesis
        from pmx import Model
        pdb = ()

