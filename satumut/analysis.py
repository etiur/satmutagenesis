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
    parser.add_argument("--thres", required=False, default=0.0, type=float,
                        help="The threshold for the improvement which will affect what will be included in the summary")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-x", "--xtc", required=False, action="store_true", help="Change the pdb format to xtc")

    args = parser.parse_args()

    return [args.inp, args.dpi, args.traj, args.out, args.plot, args.analyse,  args.cpus, args.thres,
            args.catalytic_distance, args.xtc]


class SimulationData:
    """
    A class to store data from the simulations
    """
    def __init__(self, folder, pdb=10, catalytic_dist=3.5):
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
        self.len_diff = None
        self.len = None
        self.name = basename(self.folder)

    def filtering(self):
        """
        Constructs a dataframe from all the reports in a PELE simulation folder
        """
        pd.options.mode.chained_assignment = None
        reports = []
        for files in glob("{}/report_*".format(self.folder)):
            rep = basename(files).split("_")[1]
            data = pd.read_csv(files, sep="    ", engine="python")
            data['#Task'].replace({1: rep}, inplace=True)
            data.rename(columns={'#Task': "ID"}, inplace=True)
            reports.append(data)
        self.dataframe = pd.concat(reports)
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - 99]

        # for the PELE profiles
        self.profile = self.dataframe.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        self.profile["Type"] = [self.name for _ in range(len(self.profile.index))]
        trajectory = self.dataframe.sort_values(by="distance0.5")
        trajectory.reset_index(drop=True, inplace=True)
        trajectory.drop(["Step", 'sasaLig', 'currentEnergy'], axis=1, inplace=True)
        self.trajectory = trajectory.iloc[:self.pdb]
        frequency = trajectory.loc[trajectory["distance0.5"] <= self.catalytic]  # frequency of catalytic poses
        self.frequency = frequency["distance0.5"].copy()
        self.len = len(frequency)
        self.binding = frequency["Binding Energy"].copy()
        self.binding.sort_values(inplace=True)
        self.binding.reset_index(drop=True, inplace=True)
        if "original" in self.folder:
            self.dist_ori = self.frequency.median()
            self.bind_ori = self.binding.median()

    def set_distance(self, ori_distance):
        """
        Set the distance difference with the mean

        Parameters
        __________
        original_distance: int
            The distance for the wild type
        """
        self.dist_diff = self.frequency - ori_distance  # improvement over the wild type catalytic distance

    def set_binding(self, ori_binding):
        """
        Set the binding energy difference

        Parameters
        __________
        original_binding: int
            The binding energy for the wild type
        """
        self.bind_diff = self.binding - ori_binding


def analyse_all(folders, wild, res_dir, position_num, traj=10, cata_dist=3.5):
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

    Returns
    ----------
    data_dict: dict
        Dictionary of SimulationData objects
    """
    data_dict = {}
    len_dict = {}
    original = SimulationData(wild, pdb=traj, catalytic_dist=cata_dist)
    original.filtering()
    data_dict["original"] = original
    len_dict["original"] = original.len
    for folder in folders:
        name = basename(folder)
        data = SimulationData(folder, pdb=traj, catalytic_dist=cata_dist)
        data.filtering()
        data.set_distance(original.dist_ori)
        data.set_binding(original.bind_ori)
        data_dict[name] = data
        len_dict[name] = data.len
    frame = pd.DataFrame(pd.Series(len_dict), columns=["frequency"])
    try:
        frame["ratio"] = frame["frequency"] / frame["frequency"].loc["original"]
    except ZeroDivisionError:
        pass
    if not os.path.exists("{}_results".format(res_dir)):
        os.makedirs("{}_results".format(res_dir))
    frame.to_csv("{}_results/freq_{}.csv".format(res_dir, position_num))
    return data_dict


def box_plot(res_dir, data_dict, position_num, dpi=800, cata_dist=3.5):
    """
    Creates a box plot of the 19 mutations from the same position

    Parameters
    ___________
    res_dir: str
        name of the results folder
    data_dict: dict
        A dictionary that contains SimulationData objects from the simulation folders
    position_num: str
        Position at the which the mutations occurred
    dpi: int, optional
        The quality of the plots produced
    """
    if not os.path.exists("{}_results/Plots/box".format(res_dir)):
        os.makedirs("{}_results/Plots/box".format(res_dir))
    # create a dataframe with only the distance differences for each simulation
    plot_dict_bind = {}
    plot_dict_freq = {}
    plot_dif_dist = {}
    plot_dif_bind = {}

    for key, value in data_dict.items():
        plot_dict_bind[key] = value.binding
        plot_dict_freq[key] = value.frequency
        if "original" not in key:
            plot_dif_bind[key] = value.bind_diff
            plot_dif_dist[key] = value.dist_diff

    dif_dist = pd.DataFrame(plot_dif_dist)
    dif_bind = pd.DataFrame(plot_dif_bind)
    data_bind = pd.DataFrame(plot_dict_bind)
    data_freq = pd.DataFrame(plot_dict_freq)

    # Distance difference boxplot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(data=dif_dist, kind="violin", palette="Accent", height=4.5, aspect=2.3, inner="quartile")
    ax.set(title="{} distance variation with respect to wild type".format(position_num))
    ax.set_ylabels("Distance variation", fontsize=8)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("{}_results/Plots/box/{}_distance_dif.png".format(res_dir, position_num), dpi=dpi)

    # Binding energy difference Box plot
    ex = sns.catplot(data=dif_bind, kind="violin", palette="Accent", height=4.5, aspect=2.3, inner="quartile")
    ex.set(title="{} binding energy variation with respect to wild type".format(position_num))
    ex.set_ylabels("Binding energy variation", fontsize=8)
    ex.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ex.set_xticklabels(fontsize=6)
    ex.set_yticklabels(fontsize=6)
    ex.savefig("{}_results/Plots/box/{}_binding_dif.png".format(res_dir, position_num), dpi=dpi)
    plt.close("all")

    # frequency boxplot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(data=data_freq, kind="violin", palette="Accent", height=4.5, aspect=2.3, inner="quartile")
    ax.set(title="{} distances less than {}".format(position_num, cata_dist))
    ax.set_ylabels("Distances", fontsize=8)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("{}_results/Plots/box/{}_distance.png".format(res_dir, position_num), dpi=dpi)

    # Binding energy boxplot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(data=data_bind, kind="violin", palette="Accent", height=4.5, aspect=2.3, inner="quartile")
    ax.set(title="{} binding energy ".format(position_num))
    ax.set_ylabels("Binding energy", fontsize=8)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("{}_results/Plots/box/{}_binding.png".format(res_dir, position_num), dpi=dpi)


def pele_profile_single(key, mutation, res_dir, wild, type_, position_num, dpi=800, mode="results"):
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
    ax = sns.relplot(x=type_, y='Binding Energy', hue="Type", style="Type", palette="muted", data=cat,
                     height=3.5, aspect=1.5, s=20, linewidth=0, alpha=0.2)

    ax.set(title="{} scatter plot of binding energy vs {} ".format(key, type_))
    ax.savefig("{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, type_,
                                                                 key, type_), dpi=dpi)
    plt.close(ax.fig)


def pele_profiles(type_, res_dir, data_dict, position_num, dpi=800, mode="results"):
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
    """
    for key, value in data_dict.items():
        if "original" not in key:
            pele_profile_single(key, value, res_dir=res_dir, wild=data_dict["original"],
                                type_=type_, position_num=position_num, dpi=dpi, mode=mode)


def all_profiles(res_dir, data_dict, position_num, dpi=800, mode="results"):
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
    """
    types = ["distance0.5", "sasaLig", "currentEnergy"]
    for type_ in types:
        pele_profiles(type_, res_dir, data_dict, position_num, dpi, mode=mode)


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


def create_report(res_dir, mutation, position_num, output="summary", analysis="distance", cata_dist=3.5, mode="results"):
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
        freq = val.len
        message = 'Mutation {}: median distance increment {}, median binding energy increment {}'.format(key, dis, bind)
        message2 = "{} accepted steps with a distance less than {} angstroms" .format(freq, cata_dist)
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
        plot3 = "{}_{}/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, mode, position_num, "currentEnergy", mut,
                                                             "currentEnergy")
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
                       cata_dist=3.5, mode="results"):
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
            "{} mutations at position {} decrease {} by {} or less".format(count, position_num, analysis, thres))
        create_report(res_dir, mutation_dict, position_num, output, analysis, cata_dist, mode=mode)
    else:
        log.warning("No mutations at position {} decrease {} by {} or less".format(position_num, analysis, thres))


def consecutive_analysis(file_name, wild=None, dpi=800, traj=10, output="summary", plot_dir=None, opt="distance",
                         cpus=10, thres=0.0, cata_dist=3.5, xtc=False):
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
        plot_dir = basename(dirname(dirname(plot_dir))).replace("_mut", "")
    for folders in pele_folders:
        base = basename(folders[0])[:-1]
        data_dict = analyse_all(folders, wild, plot_dir, base, traj=traj, cata_dist=cata_dist)
        box_plot(plot_dir, data_dict, base, dpi, cata_dist)
        all_profiles(plot_dir, data_dict, base, dpi)
        extract_all(plot_dir, data_dict, folders, cpus=cpus, xtc=xtc)
        find_top_mutations(plot_dir, data_dict, base, output, analysis=opt, thres=thres, cata_dist=cata_dist)


def main():
    inp, dpi, traj, out, folder, analysis, cpus, thres, cata_dist, xtc = parse_args()
    consecutive_analysis(inp, dpi=dpi, traj=traj, output=out, plot_dir=folder, opt=analysis, cpus=cpus, thres=thres,
                         cata_dist=cata_dist, xtc=xtc)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
