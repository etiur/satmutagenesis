"""
This script is used to analyse the results of the simulations
"""

from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename, dirname, abspath, commonprefix
import os
import sys
import re
from fpdf import FPDF
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
from .helper import Log, commonlist, find_log, isiterable
import mdtraj as md
import numpy as np
from collections import namedtuple
plt.switch_backend('agg')


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--inp", required=True,
                        help="Include the path of the completed_simulations.log created after satumut")
    parser.add_argument("--dpi", required=False, default=800, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("--traj", required=False, default=5, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--out", required=False, default="summary",
                        help="Name of the summary file created at the end of the analysis")
    parser.add_argument("--plot", required=False, help="Path of the plots folder")
    parser.add_argument("--analyse", required=False, choices=("energy", "distance", "both"), default="distance",
                        help="The metric to measure the improvement of the system")
    parser.add_argument("--cpus", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--thres", required=False, default=0.0, type=float,
                        help="The threshold for the improvement which will affect what will be included in the summary")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-x", "--xtc", required=False, action="store_true", help="Change the pdb format to xtc")
    parser.add_argument("-ex", "--extract", required=False, type=int, help="The number of steps to analyse")
    parser.add_argument("-en", "--energy_threshold", required=False, type=int,
                        help="An energy threshold that limits the points of scatter plots")
    parser.add_argument("-pw", "--profile_with", required=False, choices=("Binding Energy", "currentEnergy"),
                        default="Binding Energy", help="The metric to generate the pele profiles with")
    parser.add_argument("-at", "--atoms", required=False, nargs="+",
                        help="Series of atoms of the residues to follow in this format -> chain ID:position:atom name")
    parser.add_argument("-w", "--wild", required=False, default=None,
                        help="The path to the folder where the reports from wild type simulation are")
    args = parser.parse_args()

    return [args.inp, args.dpi, args.traj, args.out, args.plot, args.analyse,  args.cpus, args.thres,
            args.catalytic_distance, args.xtc, args.extract, args.energy_threshold, args.profile_with, args.atoms,
            args.wild]


class SimulationData:
    """
    A class to store data from the simulations in dictionaries
    """
    def __init__(self, folder, pdb=5, catalytic_dist=3.5, energy_thres=None, extract=None):
        """
        Initialize the SimulationData Object

        Parameters
        ___________
        folder: str
            path to the simulation folder
        points: int, optional
            Number of points to consider for the barplots
        pdb: int, optional
            how many pdbs to extract from the trajectories
        energy_thres: int, optional
            The binding energy to consider for catalytic poses
        data: pd.Dataframe
            A dataframe object containing the information from the reports
        extract: int, optional
            The number of steps to analyse
        """
        self.folder = folder
        self.extract = extract
        self.dataframe = None
        self.dist_diff = None
        self.profile = None
        self.trajectory = None
        self.pdb = pdb
        self.binding = None
        self.bind_diff = None
        self.catalytic = catalytic_dist
        self.frequency = None
        self.len_ratio = None
        self.len = None
        self.name = basename(folder)
        self.energy = energy_thres
        self.residence = None
        self.weight_dist = None
        self.weight_bind = None
        self.all = None
        self.followed_distance = "distance0.5"

    def filtering(self, followed_distance=None):
        """
        Generates the different
        """
        pd.options.mode.chained_assignment = None
        reports = []
        for files in glob("{}/report_*".format(self.folder)):
            residence_time = [0]
            rep = basename(files).split("_")[1]
            data = pd.read_csv(files, sep="    ", engine="python")
            data['#Task'].replace({1: rep}, inplace=True)
            data.rename(columns={'#Task': "ID"}, inplace=True)
            for x in range(1, len(data)):
                residence_time.append(data["Step"].iloc[x] - data["Step"].iloc[x-1])
            data["residence time"] = residence_time
            reports.append(data)
        self.dataframe = pd.concat(reports)
        if self.extract:
            self.dataframe = self.dataframe[self.dataframe["Step"] <= self.extract]
        self.dataframe.sort_values(by="currentEnergy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - min(int(len(self.dataframe)*0.1), 20)]
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - int(len(self.dataframe)*0.2)] # eliminating the 20% with the highest biding energies
        if followed_distance:
            self.followed_distance = followed_distance
        # extracting trajectories
        trajectory = self.dataframe.sort_values(by=self.followed_distance)
        trajectory.reset_index(drop=True, inplace=True)
        self.trajectory = trajectory.iloc[:self.pdb]
        if not self.energy:
            frequency = trajectory.loc[trajectory[self.followed_distance] <= self.catalytic]  # frequency of catalytic poses
            self.profile = self.dataframe.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        else:
            frequency = trajectory.loc[(trajectory[self.followed_distance] <= self.catalytic) &
                                       (trajectory["Binding Energy"] <= self.energy)]
            self.profile = frequency.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        # for the PELE profiles
        self.profile["Type"] = [self.name for _ in range(len(self.profile.index))]
        frequency["Type"] = [self.name for _ in range(len(frequency.index))]
        # binning
        self.all = pd.DataFrame(np.repeat(frequency[[self.followed_distance, "Binding Energy", "residence time", "Type"]].values,
                                          frequency["residence time"].values, axis=0),
                                columns=[self.followed_distance, "Binding Energy", "residence time", "Type"])


def bar_plot(res_dir, position_num, bins, interval, dpi=800, bin_type="distance", follow="distance0.5"):
    """
    Creates a box plot of the 19 mutations from the same position

    Parameters
    ___________
    res_dir: str
        name of the results folder
    position_num: str
        Position at the which the mutations occurred
    bins: tuple(pd.Dataframe, pd.Dataframe)
        A tuple of the dataframes from which the bar plots will be generated
    interval: str
        The interval used to generate the bins
    dpi: int, optional
        The quality of the plots produced
    bin_type: str
        2 possible types, distance bins and energy bins
    """
    if not os.path.exists("{}_results/Plots/bar".format(res_dir)):
        os.makedirs("{}_results/Plots/bar".format(res_dir))
    # create bar plots with each of the mutants
    median_bin, len_bin = bins

    # median bar plot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    median_bin.reset_index(inplace=True)
    median_bin.plot(x="index", kind="bar", stacked=False)
    if "distance" in bin_type:
        plt.xlabel('Energy intervals')
        plt.title("Median bar plot of {} bin -interval {}- with varying energy".format(follow, interval))
    else:
        plt.xlabel('Distance intervals')
        plt.title("Median bar plot of {} bin -interval {}- with varying {}".format(bin_type, interval, follow))

    plt.legend(loc='best')
    plt.xticks(rotation=40, fontsize=8)
    plt.tight_layout()
    plt.savefig("{}_results/Plots/bar/{}_median_{}_{}.png".format(res_dir, position_num, bin_type, follow), dpi=dpi)
    plt.close()

    # len bar plot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    len_bin.reset_index(inplace=True)
    len_bin.plot(x="index", kind="bar", stacked=False)
    if "distance" in bin_type:
        plt.xlabel('Energy intervals')
        plt.title("Frequency bar plot of {} -interval {}- with varying energy".format(bin_type, interval))
    else:
        plt.xlabel('Distance intervals')
        plt.title("Frequency bar plot of {} -interval {}- with varying {}".format(bin_type, interval, follow))

    plt.legend(loc='best')
    plt.xticks(rotation=40, fontsize=8)
    plt.tight_layout()
    plt.savefig("{}_results/Plots/bar/{}_frequency_{}_{}.png".format(res_dir, position_num, bin_type, follow), dpi=dpi)
    plt.close()


def binning(data_dict, res_dir, position_number, dpi=800, follow="distance0.5"):
    """
    Bins the values as to have a better analysis of the pele reports

    Parameters
    ___________
    data_dict: dict
        A dictionary containing the mutations as keys and the SimulationData dataframe as values
    res_dir: str
        The directory for the results
    position_num: str
        The position at the which the mutations was produced
    dpi: int
        The quality of the plots produced
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    bin_dict = {}
    for key, value in data_dict.items():
        bin_dict[key] = value.all[["Binding Energy", follow, "Type"]].copy()
    data = pd.concat(bin_dict.values())
    tup = namedtuple("bins", ["e_median", "e_len", "e_interval", "d_median", "d_len", "d_interval"])
    # creating the intervals or bins
    energy_bin = np.linspace(min(data["Binding Energy"]), min(max(data["Binding Energy"]), min(data["Binding Energy"])+3), num=5)
    distance_bin = np.linspace(min(data[follow]), min(max(data[follow]), min(data[follow])+3), num=5)
    energybin_labels = ["({}, {}]".format(round(energy_bin[i], 2), round(energy_bin[i + 1], 2)) for i in range(len(energy_bin) - 1)]
    distancebin_labels = ["({}, {}]".format(round(distance_bin[i], 2), round(distance_bin[i + 1], 2)) for i in range(len(distance_bin) - 1)]
    # The best distance with different energies
    best_distance = [data[(data["Binding Energy"].apply(lambda x: x in pd.Interval(energy_bin[i], energy_bin[i+1]))) &
                     (data[follow].apply(lambda x: x in pd.Interval(distance_bin[0], distance_bin[1])))] for i in range(len(energy_bin)-1)]
    # The best energies with different distances
    best_energy = [data[(data["Binding Energy"].apply(lambda x: x in pd.Interval(energy_bin[0], energy_bin[1]))) &
                   (data[follow].apply(lambda x: x in pd.Interval(distance_bin[i], distance_bin[i+1])))] for i in range(len(distance_bin)-1)]
    # For each bin in best_distance, I calculate the frequency, the median distance and energy for each of the mutations
    distance_len = [{key: len(frame[frame["Type"] == key]) for key in bin_dict.keys()} for frame in best_distance]
    distance_median = [{key: frame[frame["Type"] == key][follow].median() for key in bin_dict.keys()} for frame in best_distance]
    distance_energy = [{key: frame[frame["Type"] == key]["Binding Energy"].median() for key in bin_dict.keys()} for frame in best_distance]
    # For each bin in best_energy, I calculate the frequency, the median energy and distance for each of the mutations
    energy_len = [{key: len(frame[frame["Type"] == key]) for key in bin_dict.keys()} for frame in best_energy]
    energy_median = [{key: frame[frame["Type"] == key]["Binding Energy"].median() for key in bin_dict.keys()} for frame in best_energy]
    energy_distance = [{key: frame[frame["Type"] == key][follow].median() for key in bin_dict.keys()} for frame in best_energy]
    # For the energy bins, distance changes so using distance labels
    energy_median = pd.DataFrame(energy_median, index=["{} energy median".format(x) for x in distancebin_labels])
    energy_len = pd.DataFrame(energy_len, index=["{} energy freq".format(x) for x in distancebin_labels])
    energy_median = energy_median.fillna(0)
    energy_distance = pd.DataFrame(energy_distance, index=["{} energy distance median".format(x) for x in distancebin_labels])
    energy_distance = energy_distance.fillna(0)
    # For the distance bins, energy changes so using energy labels
    distance_median = pd.DataFrame(distance_median, index=["{} distance median".format(x) for x in energybin_labels])
    distance_median = distance_median.fillna(0)
    distance_len = pd.DataFrame(distance_len, index=["{} distance freq".format(x) for x in energybin_labels])
    distance_energy = pd.DataFrame(distance_energy, index=["{} distance energy median".format(x) for x in energybin_labels])
    distance_energy = distance_energy.fillna(0)

    # plotting
    bar_plot(res_dir, position_number, (distance_median.copy(), distance_len.copy()), distancebin_labels[0], dpi, "distance", follow)
    bar_plot(res_dir, position_number, (energy_median.copy(), energy_len.copy()), energybin_labels[0], dpi, "energy", follow)

    # concatenate everything
    everything = pd.concat([energy_median, energy_distance, energy_len, distance_median, distance_energy, distance_len])

    # To csv
    if not os.path.exists("{}_results/csv".format(res_dir)):
        os.makedirs("{}_results/csv".format(res_dir))
    everything.to_csv("{}_results/csv/binning_{}_{}.csv".format(res_dir, position_number, follow))
    return tup(energy_median, energy_len, energybin_labels, distance_median, distance_len, distancebin_labels)


def analyse_all(folders, wild, traj=5, cata_dist=3.5, energy_thres=None, extract=None, follow="distance0.5"):
    """
    Analyse all the 19 simulations folders and build SimulationData objects for each of them

    Parameters
    ___________
    folders: list[str]
        List of paths to the different reports to be analyzed
    wild: str
        Path to the simulations of the wild type
    position_num: str
        Position at the which the mutations occurred
    traj: int, optional
        How many snapshots to extract from the trajectories
    cata_dist: float, optional
        The catalytic distance
    extract: int, optional
        The number of steps to analyse
    energy_thres: int, optional
        The binding energy to consider for catalytic poses

    Returns
    ----------
    data_dict: dict
        Dictionary of SimulationData objects
    """
    data_dict = {}
    # run SimulationDAata for each of the mutant simulations
    original = SimulationData(wild, pdb=traj, catalytic_dist=cata_dist, energy_thres=energy_thres, extract=extract)
    original.filtering(follow)
    data_dict["original"] = original
    for folder in folders:
        name = basename(folder)
        data = SimulationData(folder, pdb=traj, catalytic_dist=cata_dist, energy_thres=energy_thres, extract=extract)
        data.filtering(follow)
        data_dict[name] = data

    return data_dict


def pele_profile_single(key, mutation, res_dir, wild, type_, position_num, dpi=800, mode="results",
                        profile_with="Binding Energy", follow="distance0.5"):
    """
    Creates a plot for a single mutation
    Parameters
    ___________
    key: str
        name for the axis title and plot
    mutation: SimulationData
        A SimulationData object
    res_dir: str
        name of the results folder
    wild: SimulationData
        SimulationData object that stores data for the wild type protein
    type_: str
        Type of scatter plot - distance0.5, sasaLig or currentEnergy
    position_num: str
        name for the folder to keep the images from the different mutations
    dpi: int, optional
        Quality of the plots
    profile_with: str, optional
        The metric to generate the pele profiles with
    """
    # Configuring the plot
    sns.set(font_scale=1.2)
    sns.set_style("ticks")
    sns.set_context("paper")
    original = wild.profile
    distance = mutation.profile
    cat = pd.concat([distance, original], axis=0)
    cat_1 = pd.concat([original, distance], axis=0)
    # Creating the scatter plots
    if not os.path.exists("{}_{}/Plots/{}/scatter_{}_{}".format(res_dir, mode, follow, position_num, type_)):
        os.makedirs("{}_{}/Plots/{}/scatter_{}_{}".format(res_dir, mode, follow, position_num, type_))

    ax = sns.relplot(x=type_, y=profile_with, hue="Type", style="Type", sizes=(40, 400), size="residence time",
                     palette="Set2", data=cat, linewidth=0, style_order=cat["Type"].unique(),
                     hue_order=cat["Type"].unique(), height=4, aspect=2)
    ex = sns.relplot(x=type_, y=profile_with, hue="Type", style="Type", sizes=(40, 400), size="residence time",
                     palette="Set2", data=cat_1, linewidth=0, style_order=cat["Type"].unique(),
                     hue_order=cat["Type"].unique(), height=4, aspect=2)
    ax.set(title="{} scatter plot of {} vs {} ".format(key, profile_with, type_))
    ex.set(title="{} scatter plot of {} vs {} ".format(key, profile_with, type_))
    ax.savefig(
        "{}_{}/Plots/{}/scatter_{}_{}/{}_{}_1.png".format(res_dir, mode, follow, position_num, type_, key, type_),
        dpi=dpi)
    ex.savefig(
        "{}_{}/Plots/{}/scatter_{}_{}/{}_{}_2.png".format(res_dir, mode, follow, position_num, type_, key, type_),
        dpi=dpi)
    plt.close(ax.fig)
    plt.close(ex.fig)


def pele_profiles(type_, res_dir, data_dict, position_num, dpi=800, mode="results", profile_with="Binding Energy",
                  follow="distance0.5"):
    """
    Creates a scatter plot for each of the 19 mutations from the same position by comparing it to the wild type

    Parameters
    ___________
    type_: str
        distance0.5, sasaLig or currentEnergy - different possibilities for the scatter plot
    res_dir: str
        Name of the results folder
    data_dict: dict
        A dictionary that contains SimulationData objects from the 19 simulation folders
    position_num: str
        Name for the folders where you want the scatter plot go in
    dpi: int, optional
        Quality of the plots
    mode: str, optional
        The name of the results folder, if results then activity mode if RS then rs mode
    profile_with: str, optional
        The metric to generate the pele profiles with
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    for key, value in data_dict.items():
        if "original" not in key:
            pele_profile_single(key, value, res_dir=res_dir, wild=data_dict["original"],
                                type_=type_, position_num=position_num, dpi=dpi, mode=mode, profile_with=profile_with,
                                follow=follow)


def all_profiles(res_dir, data_dict, position_num, dpi=800, mode="results", profile_with="Binding Energy",
                 follow="distance0.5"):
    """
    Creates all the possible scatter plots for the same mutated position

    Parameters
    ___________
    res_dir: str
        Name of the results folder for the output
    data_dict: dict
        A dictionary that contains SimulationData objects from the simulation folders
    position_num: str
        name for the folders where you want the scatter plot go in
    dpi: int, optional
        Quality of the plots
    mode: str, optional
        The name of the results folder, if results then activity mode if RS then rs mode
    profile_with: str, optional
        The metric to generate the pele profiles with
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    if profile_with == "Binding Energy":
        types = [follow, "sasaLig", "currentEnergy"]
    else:
        types = [follow, "sasaLig", "Binding Energy"]
    for type_ in types:
        pele_profiles(type_, res_dir, data_dict, position_num, dpi, mode=mode, profile_with=profile_with, follow=follow)


def extract_snapshot_xtc(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind,
                         follow="distance0.5"):
    """
    A function that extracts pdbs from xtc files

    Parameters
    ___________
    res_dir: str
        Name of the results folder where to store the output
    simulation_folder: str
        Path to the simulation folder
    f_id: str
        trajectory file ID
    position_num: str
        The folder name for the output of this function for the different simulations
    mutation: str
        The folder name for the output of this function for one of the simulations
    step: int
        The step in the trajectory you want to keep
    dist: float
        The distance between ligand and protein (used as name for the result file - not essential)
    bind: float
        The binding energy between ligand and protein (used as name for the result file - not essential)
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """

    if not os.path.exists("{}_results/{}_{}/{}_pdbs".format(res_dir, follow, position_num, mutation)):
        os.makedirs("{}_results/{}_{}/{}_pdbs".format(res_dir, follow, position_num, mutation))

    trajectories = glob("{}/*trajectory*_{}.*".format(simulation_folder, f_id))
    topology = "{}/input/{}_processed.pdb".format(dirname(dirname(simulation_folder)), mutation)
    if len(trajectories) == 0 or not os.path.exists(topology):
        sys.exit("Trajectory_{} or topology file not found".format(f_id))

    # load the trajectory and write it to pdb
    traj = md.load_xtc(trajectories[0], topology)
    name = "traj{}_step{}_dist{}_bind{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2))
    path_ = "{}_results/{}_{}/{}_pdbs".format(res_dir, follow, position_num, mutation)
    traj[int(step)].save_pdb(os.path.join(path_, name))


def extract_snapshot_from_pdb(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind,
                              follow="distance0.5"):
    """
    Extracts PDB files from trajectories

    Parameters
    ___________
    res_dir: str
        Name of the results folder where to store the output
    simulation_folder: str
        Path to the simulation folder
    f_id: str
        trajectory file ID
    position_num: str
        The folder name for the output of this function for the different simulations
    mutation: str
        The folder name for the output of this function for one of the simulations
    step: int
        The step in the trajectory you want to keep
    dist: float
        The distance between ligand and protein (used as name for the result file - not essential)
    bind: float
        The binding energy between ligand and protein (used as name for the result file - not essential)
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    if not os.path.exists("{}_results/{}_{}/{}_pdbs".format(res_dir, follow, position_num, mutation)):
        os.makedirs("{}_results/{}_{}/{}_pdbs".format(res_dir, follow, position_num, mutation))

    f_in = glob("{}/*trajectory*_{}.*".format(simulation_folder, f_id))
    if len(f_in) == 0:
        sys.exit("Trajectory_{} not found. Be aware that PELE trajectories must contain the label 'trajectory' in "
                 "their file name to be detected".format(f_id))
    f_in = f_in[0]
    with open(f_in, 'r') as res_dirfile:
        file_content = res_dirfile.read()
    trajectory_selected = re.search(r'MODEL\s+{}(.*?)ENDMDL'.format(int(step) + 1), file_content, re.DOTALL)

    # Output Snapshot
    traj = []
    path_ = "{}_results/{}_{}/{}_pdbs".format(res_dir, follow, position_num, mutation)
    name = "traj{}_step{}_dist{}_bind{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2))
    with open(os.path.join(path_, name), 'w') as f:
        traj.append("MODEL     {}".format(int(step) + 1))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))


def extract_10_pdb_single(info, res_dir, data_dict, xtc=False, follow="distance0.5"):
    """
    Extracts the top 10 distances for one mutation

    Parameters
    ___________
    info: iterable
       An iterable with the variables simulation_folder, position_num and mutation
    res_dir: str
       Name of the results folder
    data_dict: dict
       A dictionary that contains SimulationData objects from the simulation folders
    xtc: bool, optional
        Set to true if the pdb is in xtc format
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    simulation_folder, position_num, mutation = info
    data = data_dict[mutation]
    for ind in data.trajectory.index:
        ids = data.trajectory["ID"][ind]
        step = data.trajectory["numberOfAcceptedPeleSteps"][ind]
        dist = data.trajectory[follow][ind]
        bind = data.trajectory["Binding Energy"][ind]
        if not xtc:
            extract_snapshot_from_pdb(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind, follow)
        else:
            extract_snapshot_xtc(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind, follow)


def extract_all(res_dir, data_dict, folders, xtc=False, follow="distance0.5"):
    """
    Extracts the top 10 distances for the 19 mutations at the same position

    Parameters
    ___________
    res_dir: str
       name of the results folder
    data_dict: dict
       A dictionary that contains SimulationData objects from the 19 simulation folders
    folders: str
       Path to the folder that has all the simulations at the same position
    cpus: int, optional
       How many cpus to paralelize the function
    xtc: bool, optional
        Set to true if the pdb is in xtc format
    function: function
        a extract pdb function
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    args = []
    for pele in folders:
        name = basename(pele)
        output = name[:-1]
        args.append((pele, output, name))

    # paralelizing the function
    for arg in args:
        extract_10_pdb_single(arg, res_dir, data_dict, xtc=xtc, follow=follow)


def create_report(res_dir, mutation, position_num, output="summary", analysis="distance", cata_dist=3.5,
                  mode="results", profile_with="Binding Energy", follow="distance0.5"):
    """
    Create pdf files with the plots of chosen mutations and the path to the

    Parameters
    ___________
    res_dir: str
       Name of the results folder
    mutation: list
       A dictionary of SimulationData objects {key: SimulationData}
    position_num: str
       part of the path to the plots, the position that was mutated
    output: str, optional
       The pdf filename without the extension
    analysis: str, optional
       Type of the analysis (distance, binding or all)
    cata_dist: float, optional
        The catalytic distance
    follow: str, optional
        The column name of the different followed distances during PELE simulation

    Returns
    _______
    name: str
       The path of the pdf file
    """
    pdf = FPDF()
    pdf.set_top_margin(17.0)
    pdf.set_left_margin(15.0)
    pdf.set_right_margin(15.0)
    pdf.add_page()
    interval, results = mutation
    # Title
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, "Best mutations in terms of distance and/or binding energy", align='C', ln=1)
    pdf.set_font('Arial', '', size=10)
    for ind in results.index:
        message = 'Mutation {}: median distance {}, median binding energy {}'.format(ind, results["distance"].loc[ind], results["energy"].loc[ind])
        message1 = 'Intervals: {} for distance, {} for energy, Steps {}'.format(interval[1], interval[0], results["freq"].loc[ind])
        pdf.ln(3)  # linebreaks
        pdf.cell(0, 5, message1, ln=1)
        pdf.ln(3)
        pdf.cell(0, 5, message, ln=1)
        pdf.ln(8)  # linebreaks
        # Top poses
        pdf.set_font('Arial', 'B', size=12)
        pdf.cell(0, 10, "Path to the top poses", align='C', ln=1)
        pdf.set_font('Arial', size=10)
        pdf.ln(5)
        path = "{}_{}/{}_{}/{}_pdbs".format(res_dir, mode, follow, position_num, ind)
        pdf.cell(0, 10, "{}: {} ".format(ind, abspath(path)), ln=1)
        pdf.ln(8)
        # Scatter Plots
        pdf.ln(1000000)  # page break
        pdf.set_font('Arial', 'B', size=12)
        pdf.cell(0, 10, "Scatter plots", align='C', ln=1)
        pdf.set_font('Arial', '', size=10)
        pdf.ln(3)
        pdf.cell(0, 10, "Plots {}".format(ind), ln=1)
        pdf.ln(3)
        plot1 = "{}_{}/Plots/{}/scatter_{}_{}/{}_{}.png".format(res_dir, mode, follow, position_num, follow, ind,
                                                                follow)
        plot2 = "{}_{}/Plots/{}/scatter_{}_{}/{}_{}.png".format(res_dir, mode, follow, position_num, "sasaLig", ind,
                                                                "sasaLig")
        if profile_with == "Binding Energy":
            plot3 = "{}_{}/Plots/{}/scatter_{}_{}/{}_{}.png".format(res_dir, mode, follow, position_num,
                                                                    "currentEnergy", ind, "currentEnergy")
        else:
            plot3 = "{}_{}/Plots/{}/scatter_{}_{}/{}_{}.png".format(res_dir, mode, follow, position_num,
                                                                    "Binding Energy", ind, "Binding Energy")
        pdf.image(plot1, w=180)
        pdf.ln(3)
        pdf.image(plot2, w=180)
        pdf.ln(1000000)  # page break
        pdf.ln(3)
        pdf.image(plot3, w=180)
        pdf.ln(1000000)  # page break

    # box plots
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Bar plots of {} bins".format(analysis), align='C', ln=1)
    pdf.ln(8)
    if analysis == "distance":
        box1 = "{}_{}/Plots/bar/{}_frequency_{}_{}.png".format(res_dir, mode, position_num, "distance", follow)
        box2 = "{}_{}/Plots/bar/{}_median_{}_{}.png".format(res_dir, mode, position_num, "distance", follow)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
    elif analysis == "":
        if analysis == "distance":
            box1 = "{}_{}/Plots/bar/{}_frequency_{}_{}.png".format(res_dir, mode, position_num, "energy", follow)
            box2 = "{}_{}/Plots/bar/{}_median_{}_{}.png".format(res_dir, mode, position_num, "energy", follow)
            pdf.image(box1, w=180)
            pdf.ln(5)
            pdf.image(box2, w=180)
            pdf.ln(1000000)
    else:
        box1 = "{}_{}/Plots/bar/{}_frequency_{}.png".format(res_dir, mode, position_num, "distance", follow)
        box2 = "{}_{}/Plots/bar/{}_median_{}.png".format(res_dir, mode, position_num, "distance", follow)
        box5 = "{}_{}/Plots/bar/{}_median_{}.png".format(res_dir, mode, position_num, "energy", follow)
        box4 = "{}_{}/Plots/bar/{}_frequency_{}.png".format(res_dir, mode, position_num, "energy", follow)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
        pdf.image(box4, w=180)
        pdf.ln(5)
        pdf.image(box5, w=180)
        pdf.ln(1000000)

    # Output report
    name = "{}_{}/{}_{}_{}.pdf".format(res_dir, mode, output, position_num, follow)
    pdf.output(name, 'F')
    return name


def find_top_mutations(res_dir, bins, position_num, output="summary", analysis="distance", thres=0.0,
                       cata_dist=3.5, mode="results", energy_thres=None, profile_with="Binding Energy",
                       follow="distance0.5"):
    """
    Finds those mutations that decreases the binding distance and binding energy and creates a report

    Parameters
    ___________
    res_dir: str
       Name of the results folder
    data_dict: dict
       A dictionary of SimulationData objects that holds information for all mutations
    position_num: str
       The position that was mutated
    output: str, optional
       Name of the reports created
    analysis: str, optional
       Choose between ("distance", "binding" or "all") to specify how to filter the mutations to keep
    thres: float, optional
       Set the threshold for those mutations to be included in the pdf
    cata_dist: float, optional
        The catalytic distance
    energy_thres: int, optional
        The binding energy to consider for catalytic poses
    """
    # Find top mutations
    log = Log("{}_{}/analysis".format(res_dir, mode))
    mutation_dict = []
    # unzip the different dataframes
    e_labels = bins.e_interval
    e_median = bins.e_median
    drop_em = e_median.drop(["original"], axis=1)
    ori_em = e_median["original"]
    d_labels = bins.d_interval
    d_median = bins.d_median
    drop_dm = d_median.drop(["original"], axis=1)
    ori_dm = d_median["original"]
    d_len = bins.d_len
    # Analyse the bins
    median1 = drop_dm.loc["{} distance median".format(e_labels[0])]  # I will get for the different mutations, the distance median of the best energy and distance bin
    ene_med = drop_em.loc["{} energy median".format(d_labels[0])]  # I will get for the different mutations, the energy median of the best energy and distance bin
    if analysis == "distance":
        if ori_dm.loc["{} distance median".format(e_labels[0])] == 0:
            median1 = median1[median1 > 0]
        else:
            median1 = median1[(median1 - ori_dm.loc["{} distance median".format(e_labels[0])]) < thres]
        if not median1.empty:  # if the any of the mutations has a distance median < thres, then I concat the median distance, the frequency and the energy median
            cat = pd.concat([median1, d_len.loc["{} distance freq".format(e_labels[0])].loc[median1.index],
                             e_median.loc["{} energy median".format(d_labels[0])].loc[median1.index]], axis=1)
            cat.columns = ["distance", "freq", "energy"]
            mutation_dict = [("{} distance median".format(e_labels[0]), "{} energy median".format(d_labels[0])), cat.copy()]
    elif analysis == "energy":
        if ori_em.loc["{} energy median".format(d_labels[0])] == 0:
            ene_med = ene_med[ene_med > 0]
        else:
            ene_med = ene_med[(ene_med - ori_em.loc["{} energy median".format(d_labels[0])]) < thres]
        if not ene_med.empty:  # if the any of the mutations has an energy median < thres, then I concat the median distance, the frequency and the energy median
            cat = pd.concat([median1.loc["{} distance median".format(e_labels[0])].loc[ene_med.index],
                             d_len.loc["{} distance freq".format(e_labels[0])].loc[ene_med.index], ene_med], axis=1)
            cat.columns = ["distance", "freq", "energy"]
            mutation_dict = [("{} distance median".format(e_labels[0]), "{} energy median".format(d_labels[0])), cat.copy()]
    else:
        if ori_em.loc["{} energy median".format(d_labels[0])] == 0:
            ene_med = ene_med[ene_med > 0]
        else:
            ene_med = ene_med[(ene_med - ori_em.loc["{} energy median".format(d_labels[0])]) < thres]
        if ori_dm.loc["{} distance median".format(e_labels[0])] == 0:
            median1 = median1[median1 > 0]
        else:
            median1 = median1[(median1 - ori_dm.loc["{} distance median".format(e_labels[0])]) < thres]
        if not median1.empty and ene_med.empty:  # if the any of the mutations has an energy and distance median < thres, then I concat the median distance, the frequency and the energy median
            index = set(median1.index).intersection(ene_med.index)
            cat = pd.concat([median1.loc[index], d_len.loc["{} distance freq".format(e_labels[0])].loc[index],
                             ene_med.loc[index]], axis=1)
            cat.columns = ["distance", "freq", "energy"]
            mutation_dict = [("{} distance median".format(e_labels[0]), "{} energy median".format(d_labels[0])), cat.copy()]
    # Create a summary report with the top mutations
    if len(mutation_dict) != 0:
        log.info(
            "{} mutations at position {} decrease {} by {} or less "
            "when catalytic distance {} and binding energy {}".format(len(mutation_dict), position_num,analysis, thres,
                                                                      cata_dist, energy_thres))
        create_report(res_dir, mutation_dict, position_num, output, analysis, cata_dist, mode=mode,
                      profile_with=profile_with, follow=follow)
    else:
        log.warning("No mutations at position {} decrease {} by {} or less "
                    "when catalytic distance {} and binding energy {}".format(position_num, analysis, thres, cata_dist,
                                                                              energy_thres))


def complete_analysis(follow, folders, wild, base, dpi=800, traj=5, output="summary", plot_dir=None,
                      opt="distance", thres=0.0, cata_dist=3.5, xtc=False, extract=None, energy_thres=None,
                      profile_with="Binding Energy"):
    """
    A function that does a complete analysis of the simulation results

    Parameters
    ____________
    folders: list[str]
        List of the paths to the different simulations results of the mutants in the same position
    wild: str, optional
        The path to the wild type simulation
    base: str, optional
        The position mutated
    dpi : int, optional
       The quality of the plots
    box : int, optional
       how many points are used for the box plots
    traj : int, optional
       how many top pdbs are extracted from the trajectories
    output : str, optional
       name of the output file for the pdfs
    plot_dir : str
       Name for the results folder
    opt : str, optional
       choose if to analyse distance, energy or both
    cpus : int, optional
       How many cpus to use to extract the top pdbs
    thres : float, optional
       The threshold for the mutations to be included in the pdf
    cata_dist: float, optional
        The catalytic distance
    xtc: bool, optional
        Set to true if the pdb is in xtc format
    extract: int, optional
        The number of steps to analyse
    energy_thres: int, optional
        The binding energy to consider for catalytic poses
    profile_with: str, optional
        The metric to generate the pele profiles with
    atoms: list[str]
        Series of atoms of the residues to follow in this format -> chain ID:position:atom name, multiple of 2
    """
    data_dict = analyse_all(folders, wild, traj, cata_dist, energy_thres, extract=extract, follow=follow)
    bins = binning(data_dict, plot_dir, base, dpi=800, follow=follow)
    all_profiles(plot_dir, data_dict, base, dpi, profile_with=profile_with, follow=follow)
    extract_all(plot_dir, data_dict, folders, xtc=xtc, follow=follow)
    find_top_mutations(plot_dir, bins, base, output, analysis=opt, thres=thres, cata_dist=cata_dist,
                       energy_thres=energy_thres, profile_with=profile_with, follow=follow)
    return data_dict


def pooled_analysis(folders, wild, base, dpi=800, traj=5, output="summary", plot_dir=None, opt="distance",
                    cpus=10, thres=0.0, cata_dist=3.5, xtc=False, extract=None, energy_thres=None,
                    profile_with="Binding Energy", atoms=None):
    if atoms:
        col = ["distance{}.5".format(x) for x in range(len(atoms) // 2)]
    else:
        col = ["distance0.5"]

    p = mp.Pool(cpus)
    func = partial(complete_analysis, folders=folders, wild=wild, base=base, dpi=dpi, traj=traj, output=output,
                   plot_dir=plot_dir, opt=opt, thres=thres, cata_dist=cata_dist, xtc=xtc, extract=extract,
                   energy_thres=energy_thres, profile_with=profile_with)
    data_dicts = p.map(func, col, 1)
    p.close()
    p.join()

    # save the dataframe with the reports in csvs
    if not os.path.exists("{}_results/csv/{}".format(plot_dir, base)):
        os.makedirs("{}_results/csv/{}".format(plot_dir, base))
    for key, value in data_dicts[0].items():
        value.dataframe.to_csv("{}_results/csv/{}/{}.csv".format(plot_dir, base, key), header=True)


def consecutive_analysis(file_name, wild=None, dpi=800, traj=5, output="summary", plot_dir=None, opt="distance",
                         cpus=10, thres=0.0, cata_dist=3.5, xtc=False, extract=None, energy_thres=None,
                         profile_with="Binding Energy", atoms=None):
    """
    Analysis for the different positions

    Parameters
    ___________
    file_name : list[str], str
        An iterable that contains the path to the reports of the different simulations or the path to the directory
        where the simulations are
    wild: str, optional
        The path to the wild type simulation
    dpi : int, optional
       The quality of the plots
    box : int, optional
       how many points are used for the box plots
    traj : int, optional
       how many top pdbs are extracted from the trajectories
    output : str, optional
       name of the output file for the pdfs
    plot_dir : str
       Name for the results folder
    opt : str, optional
       choose if to analyse distance, energy or both
    cpus : int, optional
       How many cpus to use to extract the top pdbs
    thres : float, optional
       The threshold for the mutations to be included in the pdf
    cata_dist: float, optional
        The catalytic distance
    xtc: bool, optional
        Set to true if the pdb is in xtc format
    extract: int, optional
        The number of steps to analyse
    energy_thres: int, optional
        The binding energy to consider for catalytic poses
    profile_with: str, optional
        The metric to generate the pele profiles with
    atoms: list[str]
        Series of atoms of the residues to follow in this format -> chain ID:position:atom name, multiple of 2
    """
    if atoms is None:
        atoms = []
    assert len(atoms) % 2 == 0, "The number of atoms to follow should be multiple of 2"
    if isiterable(file_name):
        pele_folders = commonlist(file_name)
    elif os.path.exists("{}".format(file_name)):
        folder, wild = find_log(file_name)
        pele_folders = commonlist(folder)
    else:
        raise Exception("Pass a file with the path to the different folders")

    if not plot_dir:
        plot_dir = commonprefix(pele_folders[0])
        plot_dir = list(filter(lambda x: "_mut" in x, plot_dir.split("/")))
        plot_dir = plot_dir[0].replace("_mut", "")
    for folders in pele_folders:
        base = basename(folders[0])[:-1]
        pooled_analysis(folders, wild, base, dpi, traj, output, plot_dir, opt, cpus, thres, cata_dist, xtc, extract,
                        energy_thres, profile_with, atoms)


def main():
    inp, dpi, traj, out, folder, analysis, cpus, thres, cata_dist, xtc, extract, energy_thres, profile_with, atoms, wild \
        = parse_args()
    consecutive_analysis(inp, wild, dpi=dpi, traj=traj, output=out, plot_dir=folder, opt=analysis, cpus=cpus, thres=thres,
                         cata_dist=cata_dist, xtc=xtc, extract=extract, energy_thres=energy_thres,
                         profile_with=profile_with, atoms=atoms)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
