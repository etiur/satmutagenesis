"""
This module tests the pele files module
"""

from ..pele_files import CreateLaunchFiles, create_20sbatch
import pytest
import shutil
from os.path import basename, dirname
import os


@pytest.fixture()
def data_p():
    run = CreateLaunchFiles("data/test/PK2_F454T.pdb", "L", "ANL", "C:1:CU", "L:1:N1", 5, test=True)
    return run


class TestCreateFiles:
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

    def test_slurm(self, data_p):
        """
        A function that tests the slurm_creation function
        """
        data_p.slurm_creation("test")
        assert basename(data_p.slurm).split(".")[1] == "sh", "the file has the wrong format"
        with open(data_p.slurm, "r") as fi:
            assert fi.readline() == "#!/bin/bash\n", "slurm file not created"
        if os.path.exists(dirname(data_p.slurm)):
            shutil.rmtree(dirname(data_p.slurm))

    @pytest.mark.not_finished
    def test_nord(self, data_p):
        """
        A function that tests the slurm_nord function
        """
        data_p.slurm_nord("nord")
        assert basename(data_p.slurm).split(".")[1] == "sh", "the file has the wrong format"
        with open(data_p.slurm, "r") as fi:
            assert fi.readline() == "#!/bin/bash\n", "nord file not created"
        if os.path.exists(dirname(data_p.slurm)):
            shutil.rmtree(dirname(data_p.slurm))


def test_create_20sbatch():
    """
    A function to test the create_20sbatch function
    """
    slurm_files = create_20sbatch("L", "ANL", "C:1:CU", "L:1:N1", ["data/test/PK2_F454T.pdb"], test=True)

    with open(slurm_files[0], "r") as fi:
        assert fi.readline() == "#!/bin/bash\n", "slurm file not created"

    if os.path.exists(dirname(slurm_files[0])):
        shutil.rmtree(dirname(slurm_files[0]))