"""
This module tests the analysis module
"""

from ..analysis import SimulationData, analyse_all, box_plot, all_profiles, extract_snapshot_from_pdb
from ..analysis import create_report
import pandas as pd
import pytest
import os
import shutil


class TestSimulationData:
    """
    It is a class to test the SimulationData class
    """

    def test_filtering(self):
        """
        To test the filtering function in SimulationData
        """
        data = SimulationData("data/test/PELE/PELE_original")
        assert isinstance(data.dataframe, (pd.DataFrame, pd.Series)), "report is not read correctly"
        assert isinstance(data.profile, (pd.DataFrame, pd.Series))
        assert isinstance(data.trajectory, (pd.DataFrame, pd.Series))
        assert isinstance(data.distance, (pd.DataFrame, pd.Series))
        assert isinstance(data.binding, (pd.DataFrame, pd.Series))


@pytest.fixture()
def test_analyse_all():
    """
    Test the analyse_all function
    """
    data_dict = analyse_all("data/test/PELE/T454")
    assert type(data_dict) == dict, "data_dict not a dictionary"
    assert isinstance(data_dict["original"], (pd.DataFrame, pd.Series)), "There is no dataframe in the dictionary"

    return data_dict


def test_boxplot(test_analyse_all):
    """
    A function to test the boxplot function
    """
    box_plot("data/test/test", test_analyse_all, "test")
    assert os.path.exists("data/test/test_results/Plots/box/test_binding.png"), "the boxplot is not correct"
    if os.path.exists("data/test/test_results"):
        shutil.rmtree("data/test/test_results")


def test_pele_profiles(test_analyse_all):
    """
    Test the all_profiles function
    """
    all_profiles("data/test/test", test_analyse_all, "test")
    path = "data/test/test_results/Plots/scatter_test_distance0.5/{}_distance0.5.png"
    assert os.path.exists(path), "pele_profiles not correct"
    if os.path.exists(path):
        shutil.rmtree("data/test/test_results")


def test_extract_snapshot():
    """
    Test the extract_snapshot_from_pdb function
    """
    extract_snapshot_from_pdb("data/test/test", "data/test/PELE/PELE_original", 1, "test", "test", 10, -12, -3)
    path_ = "{}_results/distances_{}/{}_pdbs".format("data/test/test", "test", "test")
    name = "traj{}_step{}_dist{}_bind{}.pdb".format(1, 10, -12, -3)
    assert os.path.exists(os.path.join(path_, name)), "the trajectories has not been created"
    if os.path.exists(path_):
        shutil.rmtree("data/test/test_results")


def test_create_report(test_analyse_all):
    """
    To test the create report function
    """
    summary = create_report("data/test/plot", test_analyse_all, "T454")
    assert os.path.exists(summary), "the summary has not been created"
    if os.path.exists(summary):
        shutil.rmtree("data/test/test_results")


