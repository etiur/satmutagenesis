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
from helper import isiterable, Log, commonlist, find_log
import mdtraj as md
import numpy as np
plt.switch_backend('agg')


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--inp", required=True,
                        help="Include a file or list with the path to the folders with PELE simulations inside")
    parser.add_argument("--dpi", required=False, default=800, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("--traj", required=False, default=10, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--out", required=False, default="summary",
                        help="Name of the summary file created at the end of the analysis")
    parser.add_argument("--plot", required=False, help="Path of the plots folder")
    parser.add_argument("--analyse", required=False, choices=("energy", "distance", "both"), default="distance",
                        help="The metric to measure the improvement of the system")
    parser.add_argument("--cpus", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--thres", required=False, default=-0.1, type=float,
                        help="The threshold for the improvement which will affect what will be included in the summary")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.8, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-x", "--xtc", required=False, action="store_true", help="Change the pdb format to xtc")
    parser.add_argument("-ex", "--extract", required=False, type=int, help="The number of steps to analyse")
    parser.add_argument("-en", "--energy_threshold", required=False, type=int, help="The number of steps to analyse")
    parser.add_argument("-pw", "--profile_with", required=False, choices=("Binding Energy", "currentEnergy"),
                        default="Binding Energy", help="The metric to generate the pele profiles with")

    args = parser.parse_args()

    return [args.inp, args.dpi, args.traj, args.out, args.plot, args.analyse,  args.cpus, args.thres,
            args.catalytic_distance, args.xtc, args.extract, args.energy_threshold, args.profile_with]


class SimulationData:
    """
    A class to store data from the simulations
    """
    def __init__(self, folder, pdb=10, catalytic_dist=3.5, extract=None, energy_thres=None):
        """
        Initialize the SimulationData Object

        Parameters
        ___________
        folder: str
            path to the simulation folder
        points: int, optional
            Number of points to consider for the boxplots
        pdb: int, optional
            how many pdbs to extract from the trajectories
        extract: int, optional
            The number of steps to analyse
        energy_thres: int, optional
            The binding energy to consider for catalytic poses

        """
        self.folder = folder
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
        self.name = basename(self.folder)
        self.extract = extract
        self.energy = energy_thres
        self.residence = None
        self.weight_dist = None
        self.weight_bind = None
        self.all = None

    def filtering(self):
        """
        Constructs a dataframe from all the reports in a PELE simulation folder
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
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - 25]
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - int(len(self.dataframe)*0.2)] # eliminating the 20% with the highest biding energies

        # extracting trajectories
        trajectory = self.dataframe.sort_values(by="distance0.5")
        trajectory.reset_index(drop=True, inplace=True)
        self.trajectory = trajectory.iloc[:self.pdb]
        if not self.energy:
            frequency = trajectory.loc[trajectory["distance0.5"] <= self.catalytic]  # frequency of catalytic poses
        else:
            frequency = trajectory.loc[(trajectory["distance0.5"] <= self.catalytic) &
                                       (trajectory["Binding Energy"] <= self.energy)]
        # for the PELE profiles
        self.profile = frequency.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        self.profile["Type"] = [self.name for _ in range(len(self.profile.index))]
        # binning
        self.all = pd.DataFrame(np.repeat(self.profile[["distance0.5", "Binding Energy", "residence time", "Type"]].values,
                                          self.profile["residence time"].values, axis=0),
                                columns=["distance0.5", "Binding Energy", "residence time", "Type"])
        # for the csv
        self.residence = frequency["residence time"].sum()
        self.len = len(frequency)
        self.len_ratio = float(len(frequency)) / len(trajectory)
        self.frequency = self.all[["distance0.5", "residence time"]].copy()
        self.binding = self.all[["Binding Energy", "residence time"]].copy()
        self.binding.sort_values("Binding Energy", inplace=True)
        self.binding.reset_index(drop=True, inplace=True)
        self.binding = pd.DataFrame(np.repeat(self.binding.values, self.binding["residence time"].values, axis=0),
                                    columns=["Binding Energy", "residence time"])
        self.weight_dist = self.frequency["distance0.5"].median()
        self.weight_bind = self.binding["Binding Energy"].median()

    def set_distance(self, ori_distance):
        """
        Set the distance difference with the mean

        Parameters
        __________
        original_distance: int
            The distance for the wild type
        """
        self.dist_diff = self.frequency["distance0.5"] - ori_distance #improvement over the wild type catalytic distance

    def set_binding(self, ori_binding):
        """
        Set the binding energy difference

        Parameters
        __________
        original_binding: int
            The binding energy for the wild type
        """
        self.bind_diff = self.binding["Binding Energy"] - ori_binding


def bar_plot(res_dir, position_num, bins, interval, dpi=800, bin_type="Distance"):
    """
    Creates a box plot of the 19 mutations from the same position

    Parameters
    ___________
    res_dir: str
        name of the results folder
    data_dict: dict
        A dictionary that contains SimulationData objects from the simulation folders
    bins: tuple(pd.Dataframe, pd.Dataframe)
    position_num: str
        Position at the which the mutations occurred
    dpi: int, optional
        The quality of the plots produced
    """
    if not os.path.exists("{}_results/Plots/bar".format(res_dir)):
        os.makedirs("{}_results/Plots/bar".format(res_dir))
    # create bar plots with each of the mutants
    median_bin, len_bin = bins
    ind_median = np.array([x for x, _ in enumerate(median_bin.index)])
    ind_len = np.array([x for x, _ in enumerate(len_bin.index)])

    # median bar plot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    for num, key in enumerate(median_bin.columns):
        plt.bar(ind_median+(0.35*num), median_bin[key], width=0.35, label=key)
    plt.title("Median bar plot of {} bins - {}".format(bin_type, interval))
    plt.legend(loc='best')
    plt.xticks(ind_median+0.35*(len(median_bin.columns)-1)/len(median_bin.columns), median_bin.index, rotation=40, fontsize=8)
    plt.tight_layout()
    plt.savefig("{}_results/Plots/box/{}_median_{}.png".format(res_dir, position_num, bin_type), dpi=dpi)
    plt.close()

    # len bar plot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    for num, key in enumerate(len_bin.columns):
        plt.bar(ind_len+(0.35*num), len_bin[key], width=0.35, label=key)
    plt.title("frequency bar plot of {} bins- {}".format(bin_type, interval))
    plt.legend(loc='best')
    plt.xticks(ind_len+0.35*(len(median_bin.columns)-1)/len(median_bin.columns), len_bin.index, rotation=40, fontsize=8)
    plt.tight_layout()
    plt.savefig("{}_results/Plots/box/{}_frequency_{}.png".format(res_dir, position_num, bin_type), dpi=dpi)
    plt.close()


def binning(bin_dict):
    """
    Bins the values as to have a better analysis of the pele reports

    Parameters
    ___________
    bin_dict: dict
        A dictionary containing the mutations as keys and the SimulationData.all dataframe as values
    """
    data = pd.concat(bin_dict.values())
    energy_bin = np.linspace(min(data["Binding Energy"]), max(data["Binding Energy"]), num=5)
    distance_bin = np.linspace(min(data["distance0.5"]), max(data["distance0.5"]), num=5)
    energybin_labels = ["({}, {}]".format(round(energy_bin[i], 2), round(energy_bin[i + 1]), 2) for i in range(len(energy_bin) - 1)]
    distancebin_labels = ["({}, {}]".format(round(distance_bin[i], 2), round(distance_bin[i + 1]), 2) for i in range(len(distance_bin) - 1)]

    # The best distance with different energies
    distance_active = [data[(data["Binding Energy"].apply(lambda x: x in pd.Interval(energy_bin[i], energy_bin[i+1]))) &
                       (data["distance0.5"].apply(lambda x: x in pd.Interval(distance_bin[0], distance_bin[1])))] for i in range(len(energy_bin)-1)]
    # The best energies with different distances
    energy_active = [data[(data["Binding Energy"].apply(lambda x: x in pd.Interval(energy_bin[0], energy_bin[1]))) &
                     (data["distance0.5"].apply(lambda x: x in pd.Interval(distance_bin[i], distance_bin[i+1])))] for i in range(len(distance_bin)-1)]
    # For each bin in distance active, I calculate the frequency and the median of data points for each of the mutations
    distance_len = [{key: len(frame[frame["Type"] == key]) for key in bin_dict.keys()} for frame in distance_active]
    distance_median = [{key: frame[frame["Type"] == key]["distance0.5"].median() for key in bin_dict.keys()} for frame in distance_active]

    # For each bin in energy active, I calculate the frequency and the median of data points for each of the mutations
    energy_len = [{key: len(frame[frame["Type"] == key]) for key in bin_dict.keys()} for frame in energy_active]
    energy_median = [{key: frame[frame["Type"] == key]["Binding Energy"].median() for key in bin_dict.keys()} for frame in energy_active]

    # For the energy bins, distance changes so using distance labels
    energy_median = pd.DataFrame(energy_median, index=distancebin_labels)
    energy_len = pd.DataFrame(energy_len, index=distancebin_labels)
    energy_median.fillna(0, inplace=True)
    # For the distance bins, energy changes so using energy labels
    distance_median = pd.DataFrame(distance_median, index=energybin_labels)
    distance_median.fillna(0, inplace=True)
    distance_len = pd.DataFrame(distance_len, index=energybin_labels)

    # concatenate everything
    median = pd.concat([energy_median, distance_median])
    len_ = pd.concat([energy_len, distance_len])
    everything = pd.concat([len_, median])

    return everything


def analyse_all(folders, wild, res_dir, position_num, traj=10, cata_dist=3.5, extract=None, energy_thres=None):
    """
    Analyse all the 19 simulations folders and build SimulationData objects for each of them

    Parameters
    ___________
    folders: list[str]
        List of paths to the different reports to be analyzed
    wild: str
        Path to the simulations of the wild type
    res_dir: str, optional
        The directory where the results will be kept
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
    len_dict = {}
    bin_dict = {}
    weight_median = {}
    residence = {}
    len_ratio = {}
    original = SimulationData(wild, pdb=traj, catalytic_dist=cata_dist, extract=extract, energy_thres=energy_thres)
    original.filtering()
    data_dict["original"] = original
    len_dict["original"] = original.len
    bin_dict["original"] = original.all[["Binding Energy", "distance0.5"]].copy()
    weight_median["original"] = original.weight_dist
    residence["original"] = original.residence
    len_ratio["original"] = original.len_ratio
    for folder in folders:
        name = basename(folder)
        data = SimulationData(folder, pdb=traj, catalytic_dist=cata_dist, extract=extract, energy_thres=energy_thres)
        data.filtering()
        data.set_distance(original.weight_dist)
        data.set_binding(original.weight_bind)
        data_dict[name] = data
        len_dict[name] = data.len
        bin_dict[name] = data.all[["Binding Energy", "distance0.5"]].copy()
        weight_median[name] = data.weight_dist
        residence[name] = data.residence
        len_ratio[name] = data.len_ratio
    # different metrics
    if not os.path.exists("{}_results".format(res_dir)):
        os.makedirs("{}_results".format(res_dir))
    everything = binning(bin_dict)
    median = pd.DataFrame(pd.Series(weight_median), columns=["weighted median distance"])
    median["dist mut-wt"] = median["weighted median distance"] - median["weighted median distance"].loc["original"]
    # median["freq catalytic poses"] = pd.Series(len_dict)
    median["ratio catalytic vs total poses"] = pd.Series(len_ratio)
    median["frequency"] = pd.Series(len_dict)
    median["residence time"] = pd.Series(residence)
    median.to_csv("{}_results/distance_{}.csv".format(res_dir, position_num))
    everything.to_csv("{}_results/binning_{}.csv".format(res_dir, position_num))
    return data_dict


def pele_profile_single(key, mutation, res_dir, wild, type_, position_num, dpi=800, mode="results",
                        profile_with="Binding Energy"):
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
    cat = pd.concat([original, distance], axis=0)
    # Creating the scatter plots
    if not os.path.exists("{}_{}/Plots/scatter_{}_{}".format(res_dir, mode, position_num, type_)):
        os.makedirs("{}_{}/Plots/scatter_{}_{}".format(res_dir, mode, position_num, type_))
    ax = sns.relplot(x=type_, y=profile_with, hue="Type", style="Type", sizes=(40, 400), size="residence time",
                     palette="muted", data=cat, height=3.5, aspect=1.5, linewidth=0)

    ax.set(title="{} scatter plot of {} vs {} ".format(key, profile_with, type_))
    ax.savefig("{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, type_,
                                                                 key, type_), dpi=dpi)
    plt.close(ax.fig)


def pele_profiles(type_, res_dir, data_dict, position_num, dpi=800, mode="results", profile_with="Binding Energy"):
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
    """
    for key, value in data_dict.items():
        if "original" not in key:
            pele_profile_single(key, value, res_dir=res_dir, wild=data_dict["original"],
                                type_=type_, position_num=position_num, dpi=dpi, mode=mode, profile_with=profile_with)


def all_profiles(res_dir, data_dict, position_num, dpi=800, mode="results", profile_with="Binding Energy"):
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
    """
    if profile_with == "Binding Energy":
        types = ["distance0.5", "sasaLig", "currentEnergy"]
    else:
        types = ["distance0.5", "sasaLig", "Binding Energy"]
    for type_ in types:
        pele_profiles(type_, res_dir, data_dict, position_num, dpi, mode=mode, profile_with=profile_with)


def extract_snapshot_xtc(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind):
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
    """

    if not os.path.exists("{}_results/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)):
        os.makedirs("{}_results/distances_{}/{}_pdbs".format(res_dir, position_num, mutation))

    trajectories = glob("{}/*trajectory*_{}.*".format(simulation_folder, f_id))
    topology = "{}/input/{}_processed.pdb".format(dirname(dirname(simulation_folder)), mutation)
    if len(trajectories) == 0 or not os.path.exists(topology):
        sys.exit("Trajectory_{} or topology file not found".format(f_id))

    # load the trajectory and write it to pdb
    traj = md.load_xtc(trajectories[0], topology)
    name = "traj{}_step{}_dist{}_bind{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2))
    path_ = "{}_results/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)
    traj[int(step)].save_pdb(os.path.join(path_, name))


def extract_snapshot_from_pdb(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind):
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
    """
    if not os.path.exists("{}_results/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)):
        os.makedirs("{}_results/distances_{}/{}_pdbs".format(res_dir, position_num, mutation))

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
    path_ = "{}_results/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)
    name = "traj{}_step{}_dist{}_bind{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2))
    with open(os.path.join(path_, name), 'w') as f:
        traj.append("MODEL     {}".format(int(step) + 1))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))


def extract_10_pdb_single(info, res_dir, data_dict, xtc=False):
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
    """
    simulation_folder, position_num, mutation = info
    data = data_dict[mutation]
    for ind in data.trajectory.index:
        ids = data.trajectory["ID"][ind]
        step = data.trajectory["numberOfAcceptedPeleSteps"][ind]
        dist = data.trajectory["distance0.5"][ind]
        bind = data.trajectory["Binding Energy"][ind]
        if not xtc:
            extract_snapshot_from_pdb(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind)
        else:
            extract_snapshot_xtc(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind)


def extract_all(res_dir, data_dict, folders, cpus=10, xtc=False, function=None):
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
    """
    args = []
    for pele in folders:
        name = basename(pele)
        output = name[:-1]
        args.append((pele, output, name))

    # paralelizing the function
    if not function:
        function = extract_10_pdb_single
    p = mp.Pool(cpus)
    func = partial(function, res_dir=res_dir, data_dict=data_dict, xtc=xtc)
    p.map(func, args, 1)
    p.close()
    p.terminate()


def create_report(res_dir, mutation, position_num, output="summary", analysis="distance", cata_dist=3.5, mode="results",
                  profile_with="Binding Energy"):
    """
    Create pdf files with the plots of chosen mutations and the path to the

    Parameters
    ___________
    res_dir: str
       Name of the results folder
    mutation: dict
       A dictionary of SimulationData objects {key: SimulationData}
    position_num: str
       part of the path to the plots, the position that was mutated
    output: str, optional
       The pdf filename without the extension
    analysis: str, optional
       Type of the analysis (distance, binding or all)
    cata_dist: float, optional
        The catalytic distance

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

    # Title
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, "Best mutations in terms of distance and/or binding energy", align='C', ln=1)
    pdf.set_font('Arial', '', size=10)
    for key, val in mutation.items():
        dis = round(val.dist_diff.median(), 4)
        bind = round(val.bind_diff.median(), 4)
        freq = val.residence
        message = 'Mutation {}: median distance increment {}, median binding energy increment {}'.format(key, dis, bind)
        message2 = "{} steps with a distance less than {} angstroms" .format(freq, cata_dist)
        pdf.ln(3)  # linebreaks
        pdf.cell(0, 5, message, ln=1)
        pdf.ln(3)
        pdf.cell(0, 5, message2, ln=1)
    pdf.ln(8)  # linebreaks

    # box plots
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Box plot of {}".format(analysis), align='C', ln=1)
    pdf.ln(8)
    if analysis == "distance":
        box1 = "{}_{}/Plots/box/{}_distance_dif.png".format(res_dir, mode, position_num)
        box2 = "{}_{}/Plots/box/{}_distance.png".format(res_dir, mode, position_num)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
    elif analysis == "energy":
        box1 = "{}_{}/Plots/box/{}_binding_dif.png".format(res_dir, mode, position_num)
        box2 = "{}_{}/Plots/box/{}_binding.png".format(res_dir, mode, position_num)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
    elif analysis == "both":
        box1 = "{}_{}/Plots/box/{}_distance_dif.png".format(res_dir, mode, position_num)
        box5 = "{}_{}/Plots/box/{}_distance.png".format(res_dir, mode, position_num)
        box2 = "{}_{}/Plots/box/{}_binding_dif.png".format(res_dir, mode, position_num)
        box4 = "{}_{}/Plots/box/{}_binding.png".format(res_dir, mode, position_num)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
        pdf.image(box4, w=180)
        pdf.ln(5)
        pdf.image(box5, w=180)
        pdf.ln(1000000)

    # Plots
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Scatter plots", align='C', ln=1)
    pdf.set_font('Arial', '', size=10)
    for mut, key in mutation.items():
        pdf.ln(3)
        pdf.cell(0, 10, "Plots {}".format(mut), ln=1)
        pdf.ln(3)
        plot1 = "{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, "distance0.5", mut,
                                                             "distance0.5")
        plot2 = "{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, "sasaLig", mut, "sasaLig")
        if profile_with == "Binding Energy":
            plot3 = "{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, "currentEnergy", mut,
                                                                 "currentEnergy")
        else:
            plot3 = "{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, "Binding Energy", mut,
                                                                 "Binding Energy")
        pdf.image(plot1, w=180)
        pdf.ln(3)
        pdf.image(plot2, w=180)
        pdf.ln(1000000)  # page break
        pdf.ln(3)
        pdf.image(plot3, w=180)
        pdf.ln(1000000)  # page break

    # Top poses
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Path to the top poses", align='C', ln=1)
    pdf.set_font('Arial', size=10)
    pdf.ln(5)
    for mut, key in mutation.items():
        path = "{}_{}/distances_{}/{}_pdbs".format(res_dir, mode, position_num, mut)
        pdf.cell(0, 10, "{}: {} ".format(mut, abspath(path)), ln=1)
        pdf.ln(5)

    # Output report
    name = "{}_{}/{}_{}.pdf".format(res_dir, mode, output, position_num)
    pdf.output(name, 'F')
    return name


def find_top_mutations(res_dir, data_dict, position_num, output="summary", analysis="distance", thres=0.0,
                       cata_dist=3.5, mode="results", energy_thres=None, profile_with="Binding Energy"):
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
    count = 0
    mutation_dict = {}
    for key, value in data_dict.items():
        if "original" not in key:
            if analysis == "distance" and value.dist_diff.median() <= thres:
                mutation_dict[key] = value
                count += 1
            elif analysis == "energy" and value.bind_diff.median() <= thres:
                mutation_dict[key] = value
                count += 1
            elif analysis == "both" and value.dist_diff.median() <= thres and value.bind_diff.median() <= thres:
                mutation_dict[key] = value
                count += 1

    # Create a summary report with the top mutations
    if len(mutation_dict) != 0:
        log.info(
            "{} mutations at position {} decrease {} by {} or less "
            "when catalytic distance {} and binding energy {}".format(count, position_num,analysis, thres, cata_dist,
                                                                   energy_thres))
        create_report(res_dir, mutation_dict, position_num, output, analysis, cata_dist, mode=mode,
                      profile_with=profile_with)
    else:
        log.warning("No mutations at position {} decrease {} by {} or less "
                    "when catalytic distance {} and binding energy {}".format(position_num, analysis, thres, cata_dist,
                                                                              energy_thres))


def consecutive_analysis(file_name, wild=None, dpi=800, traj=10, output="summary", plot_dir=None, opt="distance",
                         cpus=10, thres=0.0, cata_dist=3.5, xtc=False, extract=None, energy_thres=None,
                         profile_with="Binding Energy"):
    """
    Creates all the plots for the different mutated positions

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
    """
    if isiterable(file_name):
        pele_folders = commonlist(file_name)
    elif os.path.exists("{}".format(file_name)):
        folder, wild = find_log(file_name)
        pele_folders = commonlist(folder)
    else:
        raise Exception("Pass a list of the path to the different folders")

    if not plot_dir:
        plot_dir = commonprefix(pele_folders[0])
        plot_dir = list(filter(lambda x: "_mut" in x, plot_dir.split("/")))
        plot_dir = plot_dir[0].replace("_mut", "")
    for folders in pele_folders:
        base = basename(folders[0])[:-1]
        data_dict = analyse_all(folders, wild, plot_dir, base, traj=traj, cata_dist=cata_dist, extract=extract,
                                energy_thres=energy_thres)
        box_plot(plot_dir, data_dict, base, dpi, cata_dist)
        all_profiles(plot_dir, data_dict, base, dpi, profile_with=profile_with)
        extract_all(plot_dir, data_dict, folders, cpus=cpus, xtc=xtc)
        find_top_mutations(plot_dir, data_dict, base, output, analysis=opt, thres=thres, cata_dist=cata_dist,
                           energy_thres=energy_thres, profile_with=profile_with)


def main():
    inp, dpi, traj, out, folder, analysis, cpus, thres, cata_dist, xtc, extract, energy_thres, profile_with = parse_args()
    consecutive_analysis(inp, dpi=dpi, traj=traj, output=out, plot_dir=folder, opt=analysis, cpus=cpus, thres=thres,
                         cata_dist=cata_dist, xtc=xtc, extract=extract, energy_thres=energy_thres,
                         profile_with=profile_with)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
