"""
This script is used to analyse the results of the simulations for substrate with chiral centers
"""
from glob import glob
import numpy as np
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename, dirname, commonprefix, abspath
import os
import sys
import re
import matplotlib.pyplot as plt
from helper import isiterable, commonlist, find_log, Log, map_atom_string
from analysis import extract_all, all_profiles
import mdtraj as md
from fpdf import FPDF
import Bio.PDB
from multiprocessing import Process, Queue
plt.switch_backend('agg')


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--inp", required=True,
                        help="Include a file or list with the path to the folders with PELE simulations inside")
    parser.add_argument("-ip","--initial_pdb", required=True,
                        help="Include the path of input pdb of the simulation")
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
    parser.add_argument("-da", "--dihedral_atoms", required=True, nargs="+",
                        help="The 4 atom necessary to calculate the dihedrals in format chain id:res number:atom name")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.8, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-im", "--improve", required=False, choices=("R", "S"), default="R",
                        help="The enantiomer that should improve")
    parser.add_argument("-ex", "--extract", required=False, type=int, help="The number of steps to analyse")
    parser.add_argument("-en", "--energy_threshold", required=False, type=int, help="The number of steps to analyse")
    args = parser.parse_args()

    return [args.inp, args.dpi, args.traj, args.out, args.plot, args.analyse, args.cpus, args.thres,
            args.catalytic_distance, args.xtc, args.improve, args.extract, args.dihedral_atoms, args.energy_threshold,
            args.initial_pdb]


class SimulationRS:
    """
    A class to analyse the simulation data from the enantiomer analysis
    """
    def __init__(self, folder, dihedral_atoms, input_pdb, res_dir, pdb=10, catalytic_dist=3.5, extract=None,
                 energy=None):
        """
        Initialize the SimulationRS class

        Parameters
        ----------
        folder: str
            The path to the simulation folder
        dihedral_atoms: list[str]
            The 4 atoms necessary to calculate the dihedral in the form of chain id:res number:atom name
        input_pdb: str
            Path to the initial pdb
        res_dir: str
            The directory of the results
        pdb: int, optional
            The number of pdbs to extract
        catalytic_dist: float
            The catalytic distance
        extract: int, optional
            The number of steps to extract
        energy: int
            The energy threshold to be considered catalytic
        """
        self.input_pdb = input_pdb
        self.folder = folder
        self.dataframe = None
        self.profile = None
        self.catalytic = catalytic_dist
        self.trajectory = None
        self.freq_r = None
        self.freq_s = None
        self.pdb = pdb
        self.atom = dihedral_atoms[:]
        self.binding_r = None
        self.binding_s = None
        self.distance = None
        self.bind_diff = None
        self.catalytic = catalytic_dist
        self.len = None
        self.dist_diff = None
        self.binding = None
        self.name = basename(self.folder)
        self.extract = extract
        self.topology = "{}/input/{}_processed.pdb".format(dirname(dirname(folder)), basename(folder))
        self.energy = energy
        self.res_dir = res_dir

    def _match_dist(self):
        """
        match the user coordinates to pmx PDB coordinates
        """
        for i in range(len(self.atom)):
            self.atom[i] = map_atom_string(self.atom[i], self.input_pdb, self.topology)

    def _transform_coordinates(self):
        """
        Transform the coordinates in format chain id: resnum: atom name into md.topology.select expressions
        """
        self._match_dist()
        select = []
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("topo", self.topology)
        for coord in self.atom:
            print coord
            resSeq = coord.split(":")[1]
            name = coord.split(":")[2]
            try:
                resname = structure[0][coord.split(":")[0]][int(resSeq)-1].resname
            except KeyError:
                resname = list(structure[0].get_residues())[int(resSeq)-1].resname

            select.append((resSeq, name, resname))

        return select

    def dihedral(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired dihedral metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired dihedral metric
        """
        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)
        traj = traj[int(len(traj)*0.25):]
        select = self._transform_coordinates()
        Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resn {}".format(select[0][0], select[0][1], select[0][2])))
        Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resn {}".format(select[1][0], select[1][1], select[1][2])))
        Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resn {}".format(select[2][0], select[2][1], select[2][2])))
        Atom_pair_4 = int(traj.topology.select("resSeq {} and name {} and resn {}".format(select[3][0], select[3][1], select[3][2])))
        metric_list = md.compute_dihedrals(traj, [[Atom_pair_1, Atom_pair_2, Atom_pair_3, Atom_pair_4]])
        metric_list = pd.Series(np.degrees(metric_list.flatten()))
        metric_list.to_csv("{}_RS/dihedral_angles.csv".format(self.res_dir), mode="a", header=False)

    def accelerated_dihedral(self):
        """
        Paralelizes the insert atomtype function
        """
        if not os.path.exists("{}_RS".format(self.res_dir)):
            os.makedirs("{}_RS".format(self.res_dir))
        pros = []
        traject_list = sorted(glob("{}/trajectory_*.pdb".format(self.folder)), key=lambda s: int(basename(s)[:-4].split("_")[1]))
        for traj in traject_list:
            p = Process(target=self.dihedral, args=(traj,))
            p.start()
            pros.append(p)
        for p in pros:
            p.join()

    def filtering(self):
        """
        Get all the info from the reports
        """
        pd.options.mode.chained_assignment = None
        reports = []
        # read the reports
        for files in sorted(glob("{}/report_*".format(self.folder)), key=lambda s: int(basename(s).split("_")[1])):
            residence_time = [0]
            rep = basename(files).split("_")[1]
            data = pd.read_csv(files, sep="    ", engine="python")
            data['#Task'].replace({1: rep}, inplace=True)
            data.rename(columns={'#Task': "ID"}, inplace=True)
            for x in range(1, len(data)):
                residence_time.append(data["Step"].iloc[x] - data["Step"].iloc[x-1])
            data["residence time"] = residence_time
            data = data[int(len(data)*0.25):]
            reports.append(data)
        self.dataframe = pd.concat(reports)
        # read the dihedral angles and concat with the dataframe
        self.accelerated_dihedral()
        angles = pd.read_csv("{}_RS/dihedral_angles.csv".format(self.res_dir), header=None, index_col=0)
        self.dataframe["dihedral"] = angles
        # removing unwanted values
        if self.extract:
            self.dataframe = self.dataframe[self.dataframe["Step"] <= self.extract]
        self.dataframe.sort_values(by="currentEnergy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - 25]
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - 99]
        # the frequency of steps with pro-S or pro-R configurations
        if not self.energy:
            frequency = self.dataframe.loc[self.dataframe["distance0.5"] <= self.catalytic]  # frequency of catalytic poses
        else:
            frequency = self.dataframe.loc[(self.dataframe["distance0.5"] <= self.catalytic) & (self.dataframe["Binding Energy"] <= self.energy)]
        freq_r = frequency.loc[(frequency["dihedral"] <= -40) & (frequency["distance2.5"] >= -140)]
        freq_r["Type"] = ["R" for _ in range(len(freq_r))]
        freq_s = frequency.loc[(frequency["dihedral"] >= 40) & (frequency["dihedral"] <= 140)]
        freq_s["Type"] = ["S" for _ in range(len(freq_s))]
        self.len = pd.DataFrame(pd.Series({"R": len(np.repeat(freq_r.values, freq_r["residence time"].values, axis=0)),
                                           "S": len(np.repeat(freq_s.values, freq_s["residence time"].values, axis=0))})).transpose()
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
        trajectory = frequency.sort_values(by="distance0.5")
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
        self.freq_r = freq_r[["distance0.5", "Type", "residence time"]].copy()
        self.freq_r = pd.DataFrame(np.repeat(self.freq_r.values, self.freq_r["residence time"].values, axis=0),
                                      columns=["distance0.5", "Type", "residence time"])
        self.freq_s = freq_s[["distance0.5", "Type", "residence time"]].copy()
        self.freq_r = pd.DataFrame(np.repeat(self.freq_s.values, self.freq_s["residence time"].values, axis=0),
                                   columns=["distance0.5", "Type", "residence time"])
        self.distance = pd.concat([self.freq_r, self.freq_s])
        self.distance["mut"] = ["{}".format(self.name) for _ in range(len(self.distance))]
        self.binding_r = freq_r[["Binding Energy", "Type", "residence time"]].copy()
        self.binding_r = pd.DataFrame(np.repeat(self.binding_r.values, self.binding_r["residence time"].values, axis=0),
                                      columns=["Binding Energy", "Type", "residence time"])
        self.binding_s = freq_s[["Binding Energy", "Type", "residence time"]].copy()
        self.binding_s = pd.DataFrame(np.repeat(self.binding_s.values, self.binding_s["residence time"].values, axis=0),
                                      columns=["Binding Energy", "Type", "residence time"])
        self.binding = pd.concat([self.binding_r, self.binding_s])
        self.binding["mut"] = ["{}".format(self.name) for _ in range(len(self.binding))]

        # calculate the median of the distance and energies of R and S
        self.dist_r = self.freq_r["distance0.5"].median()
        self.dist_s = self.freq_s["distance0.5"].median()
        self.bind_r = self.binding_r["Binding Energy"].median()
        self.bind_s = self.binding_s["Binding Energy"].median()
        self.median = pd.DataFrame(pd.Series({"R": self.dist_r, "S": self.dist_s})).transpose()
        self.median.index = [self.name]

    def set_distance(self, ori_dist1, ori_dist2):
        """
        Set the distance difference

        Parameters
        __________
        ori_dist1: float
            The mean distance of R of the wild type
        ori_dist2: float
            The mean distance of S of the wild type
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
        ori_bind1: float
            The mean binding energy of R of the wild type
        ori_bind2: float
            The mean binding energy of S of the wild type
        """
        bind_r = self.binding_r["Binding Energy"] - ori_bind1
        bind_r = pd.concat([bind_r, self.binding_r["Type"]], axis=1)
        bind_s = self.binding_s["Binding Energy"] - ori_bind2
        bind_s = pd.concat([bind_s, self.binding_s["Type"]], axis=1)
        self.bind_diff = pd.concat([bind_r, bind_s])
        self.bind_diff["mut"] = ["{}".format(self.name) for _ in range(len(self.bind_diff))]


def analyse_rs(folders, wild, dihedral_atoms, initial_pdb, res_dir, position_num, traj=10, cata_dist=3.5,
               improve="R", extract=None, energy=None):
    """
    Analyse all the 19 simulations folders and build SimulationData objects for each of them

    Parameters
    ----------
    folders: list[str]
        List of paths to the different reports to be analyzed
    wild: str
        Path to the simulations of the wild type
    dihedral_atoms: list[str]
        The 4 atoms of the dihedral
    initial_pdb: str
        Path to the initial pdb
    res_dir: str
        The folder where the results of the analysis will be kept
    position_num: str
        Position at the which the mutations occurred
    traj: int, optional
        How many snapshots to extract from the trajectories
    cata_dist: float, optional
        The catalytic distance
    improve: str, optional
        The enantiomer that improves
    extract: int, optional
        The number of steps to analyse
    energy: int, optional
        The energy_threshold to be considered catalytic

    Returns
    --------
    data_dict: dict
        Dictionary of SimulationData objects
    """
    choice = ["R", "S"]
    choice.remove(improve)
    data_dict = {}
    len_list = []
    median_list = []
    original = SimulationRS(wild, dihedral_atoms, initial_pdb, res_dir,
                            pdb=traj, catalytic_dist=cata_dist, extract=extract, energy=energy)
    original.filtering()
    data_dict["original"] = original
    len_list.append(original.len)
    median_list.append(original.median)
    for folder in folders:
        name = basename(folder)
        data = SimulationRS(folder, dihedral_atoms, initial_pdb, res_dir,
                            pdb=traj, catalytic_dist=cata_dist, extract=extract, energy=energy)
        data.filtering()
        data.set_distance(original.dist_r, original.dist_s)
        data.set_binding(original.bind_r, original.bind_s)
        data_dict[name] = data
        len_list.append(data.len)
        median_list.append(data.median)
    # frequency of catalytic distances
    if not os.path.exists("{}_RS".format(res_dir)):
        os.makedirs("{}_RS".format(res_dir))
    len_list = pd.concat(len_list)
    len_list["enantio excess"] = (len_list[improve] - len_list[choice[0]])/ (len_list["S"] + len_list["R"]) * 100
    # median catalytic distances
    median_list = pd.concat(median_list)
    median_list["diff_R"] = median_list["R"] - median_list["R"].loc["original"]
    median_list["diff_S"] = median_list["S"] - median_list["S"].loc["original"]
    everything = pd.concat([median_list, len_list], axis=1)
    everything.to_csv("{}_RS/dist_{}.csv".format(res_dir, position_num))

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


def create_report(res_dir, mutation, position_num, output="summary", analysis="distance", cata_dist=3.5, improve="R"):
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
        dis = round(val.dist_diff["distance0.5"][val.dist_diff["Type"] == improve].median(), 4)
        bind = round(val.bind_diff["Binding Energy"][val.dist_diff["Type"] == improve].median(), 4)
        freq_r = val.len["R"][0]
        freq_s = val.len["S"][0]
        message = 'Mutation {}: median distance increment of {} {}, median binding energy increment of {} {}'.format(key, improve, dis, improve, bind)
        message2 = "{} that are R and {} that are S with a distance less than {} angstroms" .format(freq_r, freq_s, cata_dist)
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
        box1 = "{}_RS/Plots/box/{}_distance_dif.png".format(res_dir, position_num)
        box2 = "{}_RS/Plots/box/{}_distance.png".format(res_dir, position_num)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
    elif analysis == "energy":
        box1 = "{}_RS/Plots/box/{}_binding_dif.png".format(res_dir, position_num)
        box2 = "{}_RS/Plots/box/{}_binding.png".format(res_dir, position_num)
        pdf.image(box1, w=180)
        pdf.ln(5)
        pdf.image(box2, w=180)
        pdf.ln(1000000)
    elif analysis == "both":
        box1 = "{}_RS/Plots/box/{}_distance_dif.png".format(res_dir, position_num)
        box5 = "{}_RS/Plots/box/{}_distance.png".format(res_dir, position_num)
        box2 = "{}_RS/Plots/box/{}_binding_dif.png".format(res_dir, position_num)
        box4 = "{}_RS/Plots/box/{}_binding.png".format(res_dir, position_num)
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
        plot1 = "{}_RS/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, position_num, "distance0.5", mut,
                                                             "distance0.5")
        plot2 = "{}_RS/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, position_num, "sasaLig", mut, "sasaLig")
        plot3 = "{}_RS/Plots/scatter_{}_{}/{}_{}.png".format(res_dir, position_num, "currentEnergy", mut,
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
        path = "{}_RS/distances_{}/{}_pdbs".format(res_dir, position_num, mut)
        pdf.cell(0, 10, "{}: {} ".format(mut, abspath(path)), ln=1)
        pdf.ln(5)

    # Output report
    name = "{}_RS/{}_{}.pdf".format(res_dir, output, position_num)
    pdf.output(name, 'F')
    return name


def find_top_mutations(res_dir, data_dict, position_num, output="summary", analysis="distance", thres=0.0,
                       cata_dist=3.5, improve="R", energy=None):
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
    energy: int, optional
        The energy threshold to be considered catalytic
    """
    # Find top mutations
    log = Log("{}_RS/analysis".format(res_dir))
    count = 0
    mutation_dict = {}
    for key, value in data_dict.items():
        if "original" not in key:
            if analysis == "distance" and value.dist_diff["distance0.5"][value.dist_diff["Type"] == improve].median() <= thres:
                mutation_dict[key] = value
                count += 1
            elif analysis == "energy" and value.bind_diff["Binding Energy"][value.dist_diff["Type"] == improve].median() <= thres:
                mutation_dict[key] = value
                count += 1
            elif analysis == "both" and value.dist_diff["distance0.5"][value.dist_diff["Type"] == improve].median() <= thres and value.bind_diff["Binding Energy"][value.dist_diff["Type"] == improve].median() <= thres:
                mutation_dict[key] = value
                count += 1

    # Create a summary report with the top mutations
    if len(mutation_dict) != 0:
        log.info(
            "{} mutations at position {} decrease {} {} by {} or less"
            "when catalytic distance {} and binding energy {}".format(count, position_num, improve, analysis, thres,
                                                                      cata_dist, energy))
        create_report(res_dir, mutation_dict, position_num, output, analysis, cata_dist, improve)
    else:
        log.warning("No mutations at position {} decrease {} {} by {} or less"
                    "when catalytic distance {} and binding energy {}".format(position_num, improve, analysis, thres,
                                                                              cata_dist, energy))


def consecutive_analysis_rs(file_name, dihedral_atoms, initial_pdb, wild=None, dpi=800, traj=10, output="summary",
                            plot_dir=None, opt="distance", cpus=10, thres=0.0, cata_dist=3.5, xtc=False, improve="R",
                            extract=None, energy=None):
    """
    Creates all the plots for the different mutated positions

    Parameters
    ___________
    file_name : list[str]
        An iterable that contains the path to the reports of the different simulations
    dihedral_atoms: list[str]
        The 4 atoms necessary to calculate the dihedral in the form of chain id:res number:atom name
    input_pdb: str
        Path to the initial pdb
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
    imporve: str
        The enantiomer that should improve
    extract: int, optional
        The number of steps to analyse
    energy: int, optional
        The energy_threshold to be considered catalytic
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
        data_dict = analyse_rs(folders, wild, dihedral_atoms, initial_pdb, plot_dir, base, traj=traj,
                               cata_dist=cata_dist, improve=improve, extract=extract, energy=energy)
        box_plot_rs(plot_dir, data_dict, base, dpi, cata_dist)
        all_profiles(plot_dir, data_dict, base, dpi, mode="RS")
        extract_all(plot_dir, data_dict, folders, cpus=cpus, xtc=xtc, function=extract_10_pdb_single_rs)
        find_top_mutations(plot_dir, data_dict, base, output, analysis=opt, thres=thres, cata_dist=cata_dist,
                           improve=improve, energy=energy)


def main():
    inp, dpi, traj, out, folder, analysis, cpus, thres, cata_dist, xtc, improve, extract, dihedral_atoms, energy,\
        initial_pdb= parse_args()
    consecutive_analysis_rs(inp, dihedral_atoms, initial_pdb, dpi=dpi, traj=traj, output=out, plot_dir=folder, opt=analysis,
                            cpus=cpus, thres=thres, cata_dist=cata_dist, xtc=xtc, improve=improve, extract=extract,
                            energy=energy)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
