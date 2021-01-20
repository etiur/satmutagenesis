"""
This module tests if all the modules from saturated_muatgenesis are available
"""

def test_module():
    try:
        from saturated_mutagenesis import mutate_pdb
        from saturated_muatgenesis import pele_files
        from saturated_mutagenesis import analysis
        from saturated_mutagenesis import helper

    except ImportError as e:
        raise ImportError(" the following modules are missing from saturated_mutagenesis: {}".format(e))
