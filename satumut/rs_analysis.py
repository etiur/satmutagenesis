"""
This script is used to analyse the results of the simulations
"""

from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename, dirname, commonprefix
import os
import sys
import re
import matplotlib.pyplot as plt
from helper import isiterable, commonlist, find_log
from analysis import find_top_mutations, extract_all, all_profiles
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
    parser.add_argument("--plot", required=False,
                        help="Path of the plots folder")
    parser.add_argument("--analyse", required=False, choices=("energy", "distance", "both"), default="distance",
                        help="The metric to measure the improvement of the system")
    parser.add_argument("--cpus", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--thres", required=False, default=0.0, type=float,
                        help="The threshold for the improvement which will affect what will be included in the summary")
    parser.add_argument("--r1", required=False, type=float, help="Distance for the R1")
    parser.add_argument("--r2", required=False, type=float, help="Distance for the R2")
    parser.add_argument("--s1", required=False, type=float, help="Distance for the S1")
    parser.add_argument("--s2", required=False, type=float, help="Distance for the S2")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    args = parser.parse_args()

    return [args.inp, args.dpi, args.traj, args.out, args.plot, args.analyse, args.cpus, args.thres,
            args.catalytic_distance, args.xtc, args.r1, args.r2, args.s1, args.s2]


class SimulationRS:
    """
    A class to analyse the simulation data from the enantiomer analysis
    """
    def __init__(self, folder, dist1r, dist2r, dist1s, dist2s, pdb=10, catalytic_dist=3.5):
        self.folder = folder
        self.dataframe = None
        self.profile = None
        self.catalytic = catalytic_dist
        self.trajectory = None
        self.freq_r = None
        self.freq_s = None
        self.pdb = pdb
        self.r_dist1 = dist1r
        self.r_dist2 = dist2r
        self.s_dist1 = dist1s
        self.s_dist2 = dist2s
        self.binding_r = None
        self.binding_s = None
        self.distance = None
        self.bind_diff = None
        self.catalytic = catalytic_dist
        self.len = None
        self.dist_diff = None
        self.binding = None
        self.name = basename(self.folder)

    def filtering(self):

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
        # the frequency of steps with pro-S or pro-R configurations
        frequency = self.dataframe[self.dataframe["distance0.5"] <= self.catalytic]
        freq_r = frequency.loc[(frequency["distance1.5"] <= self.r_dist1) & (frequency["distance2.5"] <= self.r_dist2)]
        freq_r["Type"] = ["R" for _ in range(len(freq_r))]
        freq_s = frequency.loc[(frequency["distance3.5"] <= self.s_dist1) & (frequency["distance4.5"] <= self.s_dist2)]
        freq_s["Type"] = ["S" for _ in range(len(freq_s))]
        self.len = pd.DataFrame(pd.Series({"R": len(freq_r.index), "S": len(freq_s.index)})).transpose()
        self.len.index = [self.name]
        # for the PELE profiles
        self.profile = self.dataframe.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        type_ = []
        for x in self.profile.index:
            if x in freq_r.index:
                type_.append("R_{}".format(self.name))
            elif x in freq_s.index:
                type_.append("S_{}".format(self.name))
            else:
                type_.append("noise")
        self.profile["Type"] = type_
        # To extract the best distances
        trajectory = self.dataframe.sort_values(by="distance0.5")
        trajectory.reset_index(inplace=True)
        trajectory.drop(["Step", 'sasaLig', 'currentEnergy'], axis=1, inplace=True)
        self.trajectory = trajectory.iloc[:self.pdb]
        orien = []
        for x in self.trajectory["index"].values:
            if x in freq_r.index:
                orien.append("R")
            elif x in freq_s.index:
                orien.append("S")
            else:
                orien.append("noise")
        self.trajectory["orientation"] = orien
        # for the violin plots
        self.freq_r = freq_r[["distance0.5", "Type"]].copy()
        self.freq_s = freq_s[["distance0.5", "Type"]].copy()
        self.distance = pd.concat([self.freq_r, self.freq_s])
        self.distance["mut"] = ["{}".format(self.name) for _ in range(len(self.distance))]
        self.binding_r = freq_r[["Binding Energy", "Type"]].copy()
        self.binding_s = freq_s[["Binding Energy", "Type"]].copy()
        self.binding = pd.concat([self.binding_r, self.binding_s])
        self.binding["mut"] = ["{}".format(self.name) for _ in range(len(self.distance))]

        if "original" in self.folder:
            self.dist_r = self.freq_r["distance0.5"].mean()
            self.dist_s = self.freq_s["distance0.5"].mean()
            self.bind_r = self.binding_r["Binding Energy"].mean()
            self.bind_s = self.binding_s["Binding Energy"].mean()

    def set_distance(self, ori_dist1, ori_dist2):
        """
        Set the distance difference

        Parameters
        __________
        original_distance: int
            The distance for the wild type
        """
        dist_r = self.freq_r["distance0.5"] - ori_dist1
        dist_r = pd.concat([dist_r, self.freq_r["Type"]], axis=1)
        dist_s = self.freq_s["distance0.5"] - ori_dist2
        dist_s = pd.concat([dist_s, self.freq_s["Type"]], axis=1)
        self.dist_diff = pd.concat([dist_r, dist_s])
        self.dist_diff["mut"] = ["{}".format(self.name) for _ in range(len(self.dist_diff))]

    def set_binding(self, ori_bind1, ori_bind2):
        """
        Set the binding energy difference

        Parameters
        __________
        original_binding: int
            The binding energy for the wild type
        """
        bind_r = self.binding_r["Binding Energy"] - ori_bind1
        bind_r = pd.concat([bind_r, self.binding_r["Type"]], axis=1)
        bind_s = self.binding_s["Binding Energy"] - ori_bind2
        bind_s = pd.concat([bind_s, self.binding_s["Type"]], axis=1)
        self.bind_diff = pd.concat([bind_r, bind_s])
        self.bind_diff["mut"] = ["{}".format(self.name) for _ in range(len(self.bind_diff))]


def analyse_rs(folders, wild, dist1r, dist2r, dist1s, dist2s, res_dir, position_num, traj=10, cata_dist=3.5):
    """
    Analyse all the 19 simulations folders and build SimulationData objects for each of them

    Parameters
    ----------
    folders: list[str]
        List of paths to the different reports to be analyzed
    wild: str
        Path to the simulations of the wild type
    dist1r: float
        The first distance for r
    dist2r: float
        The second distance for r
    dist1s: float
        The first distance for s
    dist2s: float
        The second distance for s
    res_dir: str
        The folder where the results of the analysis will be kept
    position_num: str
        Position at the which the mutations occurred
    traj: int, optional
        How many snapshots to extract from the trajectories
    cata_dist: float, optional
        The catalytic distance

    Returns
    --------
    data_dict: dict
        Dictionary of SimulationData objects
    """
    data_dict = {}
    len_list = []
    original = SimulationRS(wild, dist1r, dist2r, dist1s, dist2s, pdb=traj, catalytic_dist=cata_dist)
    original.filtering()
    data_dict["original"] = original
    len_list.append(original.len)
    for folder in folders:
        name = basename(folder)
        data = SimulationRS(folder, dist1r, dist2r, dist1s, dist2s, pdb=traj, catalytic_dist=cata_dist)
        data.filtering()
        data.set_distance(original.dist_r, original.dist_s)
        data.set_binding(original.bind_r, original.bind_s)
        data_dict[name] = data
        len_list.append(data.len)
    len_list = pd.concat(len_list)
    try:
        len_list["ratio_r"] = len_list["R"] / len_list["R"].loc["original"]
    except ZeroDivisionError:
        pass
    try:
        len_list["ratio_s"] = len_list["S"] / len_list["S"].loc["original"]
    except ZeroDivisionError:
        pass

    if not os.path.exists("{}_RS".format(res_dir)):
        os.makedirs("{}_RS".format(res_dir))
    len_list.to_csv("{}_RS/freq_{}.csv".format(res_dir, position_num))
    return data_dict


def box_plot_rs(res_dir, data_dict, position_num, dpi=800, cata_dist=3.5):
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
    if not os.path.exists("{}_RS/Plots/box".format(res_dir)):
        os.makedirs("{}_RS/Plots/box".format(res_dir))
    # create a dataframe with only the distance differences for each simulation
    plot_dict_bind = []
    plot_dict_freq = []
    plot_dif_dist = []
    plot_dif_bind = []

    for key, value in data_dict.items():
        plot_dict_bind.append(value.binding)
        plot_dict_freq.append(value.distance)
        if "original" not in key:
            plot_dif_bind.append(value.bind_diff)
            plot_dif_dist.append(value.dist_diff)

    dif_dist = pd.concat(plot_dif_dist)
    dif_bind = pd.concat(plot_dif_bind)
    data_bind = pd.concat(plot_dict_bind)
    data_freq = pd.concat(plot_dict_freq)

    # Distance difference boxplot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(x="mut", y="distance0.5", hue="Type", data=dif_dist, kind="violin", palette="Accent", split=True,
                     height=4.5, aspect=2.3, inner="quartile")
    ax.set(title="{} distance variation with respect to wild type".format(position_num))
    ax.set_ylabels("Distance variation", fontsize=8)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("{}_RS/Plots/box/{}_distance_dif.png".format(res_dir, position_num), dpi=dpi)

    # Binding energy difference Box plot
    ex = sns.catplot(x="mut", hue="Type", y="Binding Energy", data=dif_bind, kind="violin", palette="Accent",
                     height=4.5, aspect=2.3, inner="quartile", split=True)
    ex.set(title="{} binding energy variation with respect to wild type".format(position_num))
    ex.set_ylabels("Binding energy variation", fontsize=8)
    ex.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ex.set_xticklabels(fontsize=6)
    ex.set_yticklabels(fontsize=6)
    ex.savefig("{}_RS/Plots/box/{}_binding_dif.png".format(res_dir, position_num), dpi=dpi)
    plt.close("all")

    # frequency boxplot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(x="mut", hue="Type", y="distance0.5", split=True, data=data_freq, kind="violin", palette="Accent",
                     height=4.5, aspect=2.3, inner="quartile")
    ax.set(title="{} distances less than {}".format(position_num, cata_dist))
    ax.set_ylabels("Distances", fontsize=8)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("{}_RS/Plots/box/{}_distance.png".format(res_dir, position_num), dpi=dpi)

    # Binding energy boxplot
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(x="mut", hue="Type", y="Binding Energy", split=True, data=data_bind, kind="violin",
                     palette="Accent", height=4.5, aspect=2.3, inner="quartile")
    ax.set(title="{} binding energy ".format(position_num))
    ax.set_ylabels("Binding energy", fontsize=8)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=6)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("{}_RS/Plots/box/{}_binding.png".format(res_dir, position_num), dpi=dpi)


def extract_snapshot_xtc_rs(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind, orientation):
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

    if not os.path.exists("{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)):
        os.makedirs("{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mutation))

    trajectories = glob("{}/*trajectory*_{}.*".format(simulation_folder, f_id))
    topology = "{}/input/{}_processed.pdb".format(dirname(dirname(simulation_folder)), mutation)
    if len(trajectories) == 0 or not os.path.exists(topology):
        sys.exit("Trajectory_{} or topology file not found".format(f_id))

    # load the trajectory and write it to pdb
    traj = md.load_xtc(trajectories[0], topology)
    name = "traj{}_step{}_dist{}_bind{}_{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2), orientation)
    path_ = "{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)
    traj[int(step)].save_pdb(os.path.join(path_, name))


def snapshot_from_pdb_rs(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind, orientation):
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
    if not os.path.exists("{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)):
        os.makedirs("{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mutation))

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
    path_ = "{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mutation)
    name = "traj{}_step{}_dist{}_bind{}_{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2), orientation)
    with open(os.path.join(path_, name), 'w') as f:
        traj.append("MODEL     {}".format(int(step) + 1))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))


def extract_10_pdb_single_rs(info, res_dir, data_dict, xtc=False):
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
        orientation = data.trajectory["orientation"][ind]
        if not xtc:
            snapshot_from_pdb_rs(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind,
                                 orientation)
        else:
            extract_snapshot_xtc_rs(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind,
                                    orientation)


def consecutive_analysis_rs(file_name, dist1r, dist2r, dist1s, dist2s, wild=None, dpi=800, traj=10, output="summary",
                            plot_dir=None, opt="distance", cpus=10, thres=0.0, cata_dist=3.5, xtc=False):
    """
    Creates all the plots for the different mutated positions

    Parameters
    ___________
    file_name : list[str]
        An iterable that contains the path to the reports of the different simulations
    dist1r: float
        The first distance for r
    dist2r: float
        The second distance for r
    dist1s: float
        The first distance for s
    dist2s: float
        The second distance for s
    wild: str
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
    elif os.path.exists("{}_mut".format(file_name)):
        wild, folder = find_log(file_name)
        pele_folders = commonlist(file_name)
    else:
        raise Exception("Pass a list of the path to the different folders")

    if not plot_dir:
        plot_dir = commonprefix(pele_folders[0])
        plot_dir = basename(dirname(dirname(plot_dir))).replace("_mutations", "")
    for folders in pele_folders:
        base = basename(folders[0])[:-1]
        data_dict = analyse_rs(folders, wild, dist1r, dist2r, dist1s, dist2s, plot_dir, base, traj=traj,
                               cata_dist=cata_dist)
        box_plot_rs(plot_dir, data_dict, base, dpi, cata_dist)
        all_profiles(plot_dir, data_dict, base, dpi, mode="RS")
        extract_all(plot_dir, data_dict, folders, cpus=cpus, xtc=xtc, function=extract_10_pdb_single_rs)
        find_top_mutations(plot_dir, data_dict, base, output, analysis=opt, thres=thres, cata_dist=cata_dist, mode="RS")


def main():
    inp, dpi, traj, out, folder, analysis, cpus, thres, cata_dist, xtc, r1, r2, s1, s2 = parse_args()
    consecutive_analysis_rs(inp, r1, r2, s1, s2, dpi=dpi, traj=traj, output=out, plot_dir=folder, opt=analysis,
                            cpus=cpus, thres=thres, cata_dist=cata_dist, xtc=xtc)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
