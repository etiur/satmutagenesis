"""
This module tests the pele files module
"""

from ..pele_files import CreateYamlFiles, create_20sbatch
import pytest
import shutil
from os.path import basename, dirname
import os
from ..__main__ import CreateSlurmFiles


@pytest.fixture()
def data_p():
    run = CreateYamlFiles("data/test/PK2_F454T.pdb", "L", "ANL", ["C:1:CU", "L:1:N1"], 5, test=True)
    return run


@pytest.fixture()
def data_sl():
    run = CreateSlurmFiles("data/test/PK2_F454T.pdb", "L", "ANL", ["C:1:CU", "L:1:N1"], ["A:134", "A:32"], 5, test=True)
    return run


class TestCreateYamlFiles:
    """
    A class to test the CreateLaunchFiles class
    """
    def test_yaml(self, data_p):
        """
        A function to test the input_creation function
        """
        data_p.input_creation("test")
        assert basename(data_p.yaml).split(".")[1] == "yaml", "the file has the wrong format"
        with open(data_p.yaml, "r") as fi:
            assert fi.readline() == "system: '{}'\n".format(data_p.input), "yaml file not created"

        if os.path.exists(dirname(data_p.yaml)):
            shutil.rmtree(dirname(data_p.yaml))

    def test_slurm(self, data_sl):
        """
        A function that tests the slurm_creation function
        """
        data_sl.slurm_creation("test")
        assert basename(data_p.slurm).split(".")[1] == "sh", "the file has the wrong format"
        with open(data_sl.slurm, "r") as fi:
            assert fi.readline() == "#!/bin/bash\n", "slurm file not created"
        if os.path.exists(data_sl.slurm):
            os.remove(data_sl.slurm)

    @pytest.mark.not_finished
    def test_nord(self, data_sl):
        """
        A function that tests the slurm_nord function
        """
        data_sl.slurm_nord("nord")
        assert basename(data_sl.slurm).split(".")[1] == "sh", "the file has the wrong format"
        with open(data_sl.slurm, "r") as fi:
            assert fi.readline() == "#!/bin/bash\n", "nord file not created"
        if os.path.exists(data_sl.slurm):
            os.remove(data_sl.slurm)


def test_create_20sbatch():
    """
    A function to test the create_20sbatch function
    """
    yaml_files = create_20sbatch("L", "ANL", ["C:1:CU", "L:1:N1"], ["data/test/PK2_F454T.pdb"], test=True)

    with open(yaml_files[0], "r") as fi:
        assert fi.readline() == "#!/bin/bash\n", "yaml file not created"

    if os.path.exists(dirname(yaml_files[0])):
        shutil.rmtree(dirname(yaml_files[0]))