"""
This script is used to analyse the results of the simulations for substrate with chiral centers
"""
import numpy as np
import pandas as pd
import argparse
import sys
import re
import matplotlib.pyplot as plt
from .helper import match_dist, check_completed_log
from .analysis import all_profiles
import mdtraj as md
import Bio.PDB
from multiprocessing import Pool
from functools import partial
from pathlib import Path
plt.switch_backend('agg')


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("-f", "--simulation_folder", required=True,
                        help="Include the path to the folder created during the simulation")
    parser.add_argument("-ip", "--initial_pdb", required=True,
                        help="Include the path of input pdb of the simulation")
    parser.add_argument("--dpi", required=False, default=800, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("--traj", required=False, default=5, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--plot", required=False, help="Path of the plots folder")
    parser.add_argument("--cpus", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-da", "--dihedral_atoms", required=True, nargs="+",
                        help="The 4 atom necessary to calculate the dihedrals in format chain id:res number:atom name")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-ex", "--extract", required=False, type=int, help="The number of steps to analyse")
    parser.add_argument("-en", "--energy_threshold", required=False, type=int, help="The number of steps to analyse")
    parser.add_argument("-pw", "--profile_with", required=False, choices=("Binding Energy", "currentEnergy"),
                        default="Binding Energy", help="The metric to generate the pele profiles with")
    parser.add_argument("-w", "--wild", required=False, default=None,
                        help="The path to the folder where the reports from wild type simulation are")
    args = parser.parse_args()

    return [args.simulation_folder, args.dpi, args.traj, args.plot, args.cpus,
            args.catalytic_distance, args.xtc, args.extract, args.dihedral_atoms, args.energy_threshold,
            args.initial_pdb, args.profile_with, args.wild]


def dihedral(trajectory, select, topology=None):
    """
    Take the PELE simulation trajectory files and returns the list of values of the desired dihedral metric

    Parameters
    ____________
    trajectory: str
        The path to the trajectory
    select: list[tuple(position number, atom name, residue name)]
        a list of tuple containing the 3 elements for the atom seletion in mdtraj
    topology: str
        Path to the topology file if xtc format

    RETURNS
    -------
    metric_list: list
            List of values of the desired dihedral metric
    """
    if ".xtc" in trajectory:
        traj = md.load_xtc(trajectory, topology)
        num = int(Path(trajectory).name.replace(".xtc", "").split("_")[1])
    else:
        traj = md.load_pdb(trajectory)
        num = int(Path(trajectory).name.replace(".pdb", "").split("_")[1])
    traj = traj[int(len(traj)*0.10):]
    Atom_pair_1 = int(traj.topology.select(f"resSeq {select[0][0]} and name {select[0][1]} and resn {select[0][2]}"))
    Atom_pair_2 = int(traj.topology.select(f"resSeq {select[1][0]} and name {select[1][1]} and resn {select[1][2]}"))
    Atom_pair_3 = int(traj.topology.select(f"resSeq {select[2][0]} and name {select[2][1]} and resn {select[2][2]}"))
    Atom_pair_4 = int(traj.topology.select(f"resSeq {select[3][0]} and name {select[3][1]} and resn {select[3][2]}"))
    metric_list = md.compute_dihedrals(traj, [[Atom_pair_1, Atom_pair_2, Atom_pair_3, Atom_pair_4]])
    metric_list = pd.Series(np.degrees(metric_list.flatten()))
    id = pd.Series([num for _ in range(len(metric_list))])
    metric_list = pd.concat([metric_list, id], axis=1)
    metric_list.columns = ["dihedral", "id"]
    return metric_list


class SimulationRS:
    """
    A class to analyse the simulation data from the enantiomer analysis
    """
    def __init__(self, folder, dihedral_atoms, input_pdb, res_dir, pdb=5, catalytic_dist=3.5, extract=None,
                 energy=None, cpus=10):
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
        cpus: int
            The number of cpus to extract the md trajectories
        """
        self.input_pdb = input_pdb
        self.folder = Path(folder)
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
        self.name = self.folder.name
        self.extract = extract
        self.topology = self.folder.parent.parent/"input"/f"{self.folder.name}_processed.pdb"
        self.energy = energy
        self.res_dir = res_dir
        self.cpus = cpus
        self.all = None
        self.followed = "distance0.5"

    def _transform_coordinates(self):
        """
        Transform the coordinates in format chain id: resnum: atom name into md.topology.select expressions
        """
        select = []
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("topo", self.topology)
        for coord in self.atom:
            resSeq = coord.split(":")[1]
            name = coord.split(":")[2]
            try:
                resname = structure[0][coord.split(":")[0]][int(resSeq)-1].resname
            except KeyError:
                resname = list(structure[0][coord.split(":")[0]].get_residues())[0].resname
            select.append((resSeq, name, resname))

        return select

    def accelerated_dihedral(self, select):
        """
        Paralelizes the insert atomtype function
        """
        Path(f"{self.res_dir}_RS/angles").mkdir(parents=True, exist_ok=True)
        traject_list = sorted(list(self.folder.glob("trajectory_*")), key=lambda s: int(s.stem.split("_")[1]))
        # parallelize the function
        p = Pool(self.cpus)
        func = partial(dihedral, select=select, topology=self.topology)
        angles = pd.concat(list(p.map(func, traject_list)))
        angles.reset_index(drop=True, inplace=True)
        angles.to_csv(f"{self.res_dir}_RS/angles/{self.name}.csv", header=True)
        return angles

    def filtering(self, follow=None):
        """
        Get all the info from the reports
        """
        if follow:
            self.followed = follow
        pd.options.mode.chained_assignment = None
        reports = []
        select = self._transform_coordinates()
        angles = self.accelerated_dihedral(select)
        # read the reports
        for files in sorted(self.folder.glob("report_*"), key=lambda s: int(s.name.split("_")[1])):
            residence_time = [0]
            rep = int(files.name.split("_")[1])
            data = pd.read_csv(files, sep="    ", engine="python")
            data['#Task'].replace({1: rep}, inplace=True)
            data.rename(columns={'#Task': "ID"}, inplace=True)
            for x in range(1, len(data)):
                residence_time.append(data["Step"].iloc[x] - data["Step"].iloc[x-1])
            data["residence time"] = residence_time
            data = data[int(len(data)*0.10):]
            reports.append(data)
        self.dataframe = pd.concat(reports)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe["dihedral"] = angles["dihedral"]
        # removing unwanted values
        if self.extract:
            self.dataframe = self.dataframe[self.dataframe["Step"] <= self.extract]
        self.dataframe.sort_values(by="currentEnergy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - min(int(len(self.dataframe)*0.1), 20)]
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe) - 99]
        # the frequency of steps with pro-S or pro-R configurations
        if not self.energy:
            frequency = self.dataframe.loc[self.dataframe[self.followed] <= self.catalytic]  # frequency of catalytic poses
            self.profile = self.dataframe.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        else:
            frequency = self.dataframe.loc[(self.dataframe[self.followed] <= self.catalytic) & (self.dataframe["Binding Energy"] <= self.energy)]
            self.profile = frequency.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)

        freq_r = frequency.loc[(frequency["dihedral"] <= -40) & (frequency["dihedral"] >= -140)]
        freq_r["Type"] = ["R" for _ in range(len(freq_r))]
        freq_s = frequency.loc[(frequency["dihedral"] >= 40) & (frequency["dihedral"] <= 140)]
        freq_s["Type"] = ["S" for _ in range(len(freq_s))]
        self.len = pd.DataFrame(pd.Series({"R": len(np.repeat(freq_r.values, freq_r["residence time"].values, axis=0)),
                                           "S": len(np.repeat(freq_s.values, freq_s["residence time"].values, axis=0))})).transpose()
        self.len.index = [self.name]
        # for the PELE profiles
        type_ = []
        for x in self.profile["dihedral"].values:
            if -40 >= x >= -140:
                type_.append(f"R_{self.name}")
            elif 40 <= x <= 140:
                type_.append(f"S_{self.name}")
            else:
                type_.append("noise")
        self.profile["Type"] = type_

        # To extract the best distances
        trajectory = frequency.sort_values(by=self.followed)
        trajectory.reset_index(inplace=True)
        trajectory.drop(["Step", 'sasaLig', 'currentEnergy'], axis=1, inplace=True)
        self.trajectory = trajectory.iloc[:self.pdb]
        orien = []
        for x in self.trajectory["dihedral"].values:
            if -40 >= x >= -140:
                orien.append("R")
            elif 40 <= x <= 140:
                orien.append("S")
            else:
                orien.append("noise")
        self.trajectory["orientation"] = orien


def analyse_rs(folders, wild, dihedral_atoms, initial_pdb, res_dir, traj=5, cata_dist=3.5, extract=None, energy=None,
               cpus=10, follow="distance0.5"):
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
    traj: int, optional
        How many snapshots to extract from the trajectories
    cata_dist: float, optional
        The catalytic distance
    extract: int, optional
        The number of steps to analyse
    energy: int, optional
        The energy_threshold to be considered catalytic
    cpus: int, optional
        The number of processors for the md trajectories
    follow: str, optional
        The column name of the different followed distances during PELE simulation

    Returns
    --------
    data_dict: dict
        Dictionary of SimulationData objects
    """
    data_dict = {}
    atoms = match_dist(dihedral_atoms, initial_pdb, wild)
    original = SimulationRS(wild, atoms, initial_pdb, res_dir, pdb=traj, catalytic_dist=cata_dist, extract=extract,
                            energy=energy, cpus=cpus)
    original.filtering(follow)
    data_dict["original"] = original
    for folder in folders:
        name = Path(folder).name
        data = SimulationRS(folder, atoms, initial_pdb, res_dir, pdb=traj, catalytic_dist=cata_dist, extract=extract,
                            energy=energy, cpus=cpus)
        data.filtering(follow)
        data_dict[name] = data

    return data_dict


def extract_snapshot_xtc_rs(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind, orientation,
                            angle, follow="distance0.5"):
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
    angle: float
        The dihedral angle
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    Path(f"{res_dir}_RS/{follow}_{position_num}/{mutation}_pdbs").mkdir(exist_ok=True, parents=True)
    trajectories = list(Path(simulation_folder).glob(f"*trajectory*_{f_id}.*"))
    topology = Path(simulation_folder).parent.parent/"input"/f"{mutation}_processed.pdb"
    if len(trajectories) == 0 or not topology.exists():
        sys.exit(f"Trajectory_{f_id} or topology file not found")

    # load the trajectory and write it to pdb
    traj = md.load_xtc(trajectories[0], topology)
    name = f"traj{f_id}_step{step}_dist{round(dist, 2)}_bind{round(bind, 2)}_{orientation}_{round(angle, 2)}.pdb"
    path_ = Path(f"{res_dir}_RS/{follow}_{position_num}/{mutation}_pdbs")
    traj[int(step)].save_pdb(path_/name)


def snapshot_from_pdb_rs(res_dir, simulation_folder, f_id, position_num, mutation, step, dist, bind, orientation,
                         angle, follow="distance0.5"):
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
    angle: float
        The dihedral angle of the trajectory
    follow: str, optional
        The column name of the different followed distances during PELE simulation
    """
    Path(f"{res_dir}_RS/{follow}_{position_num}/{mutation}_pdbs").mkdir(exist_ok=True, parents=True)
    f_in = list(Path(simulation_folder).glob(f"*trajectory*_{f_id}.*"))
    if len(f_in) == 0:
        sys.exit(f"Trajectory_{f_id} not found. Be aware that PELE trajectories must contain the label 'trajectory' in "
                 "their file name to be detected")
    f_in = f_in[0]
    with open(f_in, 'r') as res_dirfile:
        file_content = res_dirfile.read()
    trajectory_selected = re.search(r'MODEL\s+{}(.*?)ENDMDL'.format(int(step) + 1), file_content, re.DOTALL)

    # Output Snapshot
    traj = []
    path_ = Path(f"{res_dir}_RS/{follow}_{position_num}/{mutation}_pdbs")
    name = f"traj{f_id}_step{step}_dist{round(dist, 2)}_bind{round(bind, 2)}_{orientation}_{round(angle, 2)}.pdb"
    with open(path_/name, 'w') as f:
        traj.append("MODEL     {}".format(int(step) + 1))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))


def extract_10_pdb_single_rs(info, res_dir, data_dict, xtc=False, follow="distance0.5"):
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
        dist = data.trajectory["distance0.5"][ind]
        bind = data.trajectory["Binding Energy"][ind]
        orientation = data.trajectory["orientation"][ind]
        angle = data.trajectory["dihedral"][ind]
        if not xtc:
            snapshot_from_pdb_rs(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind,
                                 orientation, angle, follow)
        else:
            extract_snapshot_xtc_rs(res_dir, simulation_folder, ids, position_num, mutation, step, dist, bind,
                                    orientation, angle, follow)


def extract_all(res_dir, data_dict, folders, cpus=10, xtc=False, follow="distance0.5"):
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
        name = Path(pele).name
        output = name[:-1]
        args.append((pele, output, name))

    # paralelizing the function
    p = Pool(cpus)
    func = partial(extract_10_pdb_single_rs, res_dir=res_dir, data_dict=data_dict, xtc=xtc, follow=follow)
    p.map(func, args, 1)
    p.close()
    p.join()


def consecutive_analysis_rs(file_name, dihedral_atoms, initial_pdb, wild=None, dpi=800, traj=5,
                            plot_dir=None, cpus=10, cata_dist=3.5, xtc=False, extract=None, energy=None,
                            profile_with="Binding Energy"):
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
    plot_dir : str
       Name for the results folder
    cpus : int, optional
       How many cpus to use to extract the top pdbs
    cata_dist: float, optional
        The catalytic distance
    xtc: bool, optional
        Set to true if the pdb is in xtc format
    extract: int, optional
        The number of steps to analyse
    energy: int, optional
        The energy_threshold to be considered catalytic
    profile_with: str, optional
        The metric to generate the pele profiles with
    """
    pele_folders, plot_dir, wild = check_completed_log(file_name, wild, plot_dir)
    for folders in pele_folders:
        base = Path(folders[0]).name[:-1]
        data_dict = analyse_rs(folders, wild, dihedral_atoms, initial_pdb, plot_dir, traj, cata_dist, extract, energy, cpus)
        all_profiles(plot_dir, data_dict, base, dpi, mode="RS", profile_with=profile_with)
        extract_all(plot_dir, data_dict, folders, cpus, xtc=xtc)


def main():
    simulation_folder, dpi, traj, folder, cpus, cata_dist, xtc,  extract, dihedral_atoms, energy, initial_pdb, \
    profile_with, wild = parse_args()
    consecutive_analysis_rs(simulation_folder, dihedral_atoms, initial_pdb, wild, dpi=dpi, traj=traj, plot_dir=folder,
                            cpus=cpus, cata_dist=cata_dist, xtc=xtc, extract=extract, energy=energy,
                            profile_with=profile_with)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
