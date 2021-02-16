"""
This module provides test of the mutate_pdb.py's functions
"""
import pytest
from ..mutate_pdb import Mutagenesis, generate_mutations
import os
from os.path import basename


@pytest.fixture()
def data_m():
    run = Mutagenesis("data/test/PK2_F454T.pdb", "A:135", "data/test")
    run._check_coords()
    return run


class TestMutagenesis:
    """
    Test the class Mutagenesis on mutate.py
    """
    def test_mutate(self, data_m):
        """
        It checks if the pmx library is working correctly
        """
        from pmx import Model

        assert isinstance(data_m.model, Model), "model not being an instance of Model in pmx"
        data_m.mutate(data_m.chain.residues[134], "SER", data_m.rotamers)
        assert data_m.chain.residues[134].resname == "SER", "The mutate function is broken"

    @pytest.mark.slow_s
    def test_insert(self, data_m):
        """
        Checks if insert_atomtype function within the Mutagenesis works
        """
        # read in user input
        file_ = data_m.single_mutagenesis("ALA")
        assert "135A" in basename(file_), "naming of the file incorrect"
        assert basename(file_).split(".")[1] == "pdb", "the file format is wrong"

        with open(data_m.input, "r") as initial:
            initial_lines = initial.readlines()
        # read in preprocessed input
        with open("data/test/{}".format(file_), "r") as prep:
            prep_lines = prep.readlines()

        for line in prep_lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                for linex in initial_lines:
                    if linex.startswith("HETATM") or linex.startswith("ATOM"):
                        assert line[66:81].strip() == linex[66:81].strip(), "atom type insertion incorrect"
                        break
                break

        if os.path.exists(file_):
            os.remove(file_)


@pytest.mark.slow
def test_generate_mutations():
    """
    Tests the generate mutations from the mutate_pdb module
    """
    pdbs = generate_mutations("data/test/PK2_F454T.pdb", ["A:135"], pdb_dir="data/test")
    assert len(pdbs) >= 19, "failure to generate all the files"
    assert basename(pdbs[0]).split(".")[1] == "pdb", "the file has a wrong format"
    for f in pdbs:
        if os.path.isfile(f):
            os.remove(f)






