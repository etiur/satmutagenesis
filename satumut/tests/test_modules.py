"""
This module tests if all the modules from saturated_muatgenesis are available
"""


def test_module():
    try:
        from .. import mutate_pdb
        from .. import pele_files
        from .. import analysis
        from .. import helper

    except ImportError as e:
        raise ImportError(" the following modules are missing from saturated_mutagenesis: {}".format(e))
