from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename
import os
import matplotlib.pyplot as plt
import sys
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--pele", required=False, default="./",
                        help="Include a file with names of the different folders with PELE simulations inside")
    parser.add_argument("--dpi", required=False, default=1000, type=int,
                        help="Set the quality of the plots")
    args = parser.parse_args()

    return args.pele, args.dpi


class SimulationData:
    def __init__(self, folder):
        """
        folder (str):  path to the simulation folder
        """
        self.folder = folder
        self.dataframe = None
        self.distance = None
        self.distribution = None
        self.profile = None
        self.trajectory = None

    def filtering(self):
        """
        Constructs a dataframe from all the reports in a PELE simulation folder with the best 20% binding energies
        and a Series with the 100 best ligand distances
        """
        reports = []
        for files in glob("{}/output/0/report_*".format(self.folder)):
            rep = basename(files).split("_")[1]
            data = pd.read_csv(files, sep="    ", engine="python")
            data['#Task'].replace({1: rep}, inplace=True)
            data.rename(columns={'#Task': "ID"}, inplace=True)
            reports.append(data)
        self.dataframe = pd.concat(reports)
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.iloc[:len(self.dataframe)-99]

        # for the PELE profiles
        self.profile = self.dataframe.drop(["Step", "numberOfAcceptedPeleSteps", 'ID'], axis=1)
        self.trajectory = self.dataframe.sort_values(by="distance0.5")
        self.trajectory.reset_index(drop=True, inplace=True)
        self.trajectory.drop(["Step", 'sasaLig', 'currentEnergy'], axis=1, inplace=True)
        self.trajectory = self.trajectory.head(10)

        # For the box plots
        data_20 = self.dataframe.head(len(self.dataframe) * 20 / 100)
        dist20 = data_20["distance0.5"].copy()
        dist20.sort_values(inplace=True)
        dist20.reset_index(drop=True, inplace=True)
        self.distance = dist20.head(min(100, len(dist20)))  # 20 - 100
        if "original" in self.folder:
            self.distance = dist20[0].copy()

    def set_distribution(self, original_distance):
        """
        Stores the difference between the mutated ligand distances and the original ligand distance
        """
        self.distribution = self.distance - original_distance


def analyse_all(folders="."):
    """
    folders (str): path to the different PELE simulation folders to be analyzed
    """
    # create a SimulationData object for each pele simulation
    data_dict = {}
    original = SimulationData("PELE_original")
    original.filtering()
    data_dict["original"] = original
    for folder in glob("{}/PELE_*".format(folders)):
        name = basename(folder)
        data = SimulationData(folder)
        data.filtering()
        data.set_distribution(original.distance)
        data_dict[name[5:]] = data

    return data_dict


def box_plot(data_dict, name, dpi=1000):
    """
    Creates a box plot of the 19 mutations from the same position
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    """
    if not os.path.exists("results/Plots/box"):
        os.makedirs("results/Plots/box")
    # create a dataframe with only the distance differences for each simulation
    plot_dict = {}
    for key, value in data_dict.items():
        if "original" not in key:
            plot_dict[key] = value.distribution
    dataframe = pd.DataFrame(plot_dict)
    # Creates a boxplot from the dataframe
    sns.set(font_scale=1.9)
    sns.set_style("ticks")
    sns.set_context("paper")
    ax = sns.catplot(data=dataframe, kind="box", palette="Accent", height=4.5, aspect=2.3)
    ax.set(title="{} distance variation with respect to wild type".format(name))
    ax.set_ylabels("Distance variation", fontsize=9)
    ax.set_xlabels("Mutations {}".format(name), fontsize=9)
    ax.set_xticklabels(fontsize=7)
    ax.set_yticklabels(fontsize=7)
    ax.savefig("results/Plots/box/{}.png".format(name), dpi=dpi)


def pele_profile_single(wild, key, types, name, mutation, dpi=1000):
    """
    Creates a plot for a single mutation
    wild (SimulationData): SimulationData object that stores data for the wild type protein
    key (str): name of the mutation
    types (str): Type of scatter plot - distance0.5, sasaLig or currentEnergy
    name (str): name for the folder to keep the images
    mutation (SimulationData): SimulationData object that stores data for the mutated protein
    """
    sns.set(font_scale=1)
    sns.set_style("ticks")
    sns.set_context("paper")
    original = wild.profile
    original.index = ["Wild type"] * len(original)
    distance = mutation.profile
    distance.index = [key] * len(distance)
    cat = pd.concat([original, distance], axis=0)
    cat.index.name = "Type"
    cat.reset_index(inplace=True)
    if types == "currentEnergy":
        if not os.path.exists("results/Plots/scatter_{}_{}/{}".format(name, types, "distance0.5")):
            os.makedirs("results/Plots/scatter_{}_{}/{}".format(name, types, "distance0.5"))
        if not os.path.exists("results/Plots/scatter_{}_{}/{}".format(name, types, "sasaLig")):
            os.makedirs("results/Plots/scatter_{}_{}/{}".format(name, types, "sasaLig"))

        norm = plt.Normalize(cat["sasaLig"].min(), cat["sasaLig"].max())
        norm2 = plt.Normalize(cat["distance0.5"].min(), cat["distance0.5"].max())
        ax = sns.relplot(x=types, y='Binding Energy', hue="sasaLig", style="Type", palette='RdBu', data=cat,
                         height=3.8, aspect=1.8, hue_norm=norm, s=100, linewidth=0)
        ex = sns.relplot(x=types, y='Binding Energy', hue="distance0.5", style="Type", palette='RdBu', data=cat,
                         height=3.8, aspect=1.8, hue_norm=norm2, s=100, linewidth=0)
        ex.set(title="{} scatter plot of binding energy vs {} ".format(key, types))
        ex.savefig("results/Plots/scatter_{}_{}/{}/{}_{}.png".format(name, types, "distance0.5", key, types), dpi=dpi)
        ax.savefig("results/Plots/scatter_{}_{}/{}/{}_{}.png".format(name, types, "sasaLig", key, types), dpi=dpi)

    else:
        if not os.path.exists("results/Plots/scatter_{}_{}".format(name, types)):
            os.makedirs("results/Plots/scatter_{}_{}".format(name, types))
        ax = sns.relplot(x=types, y='Binding Energy', hue="Type", style="Type", palette="Set1", data=cat,
                         height=3.8, aspect=1.8, s=100, linewidth=0)
        ax.set(title="{} scatter plot of binding energy vs {} ".format(key, types))
        ax.savefig("results/Plots/scatter_{}_{}/{}_{}.png".format(name, types, key, types), dpi=dpi)


def pele_profiles(data_dict, name, types, dpi=1000):
    """
    Creates a scatter plot for each of the 19 mutations from the same position by comparing it to the wild type
    data_dict (dict): A dictionary that contains SimulationData objects from the 19 simulation folders
    name (str): name for the folders where you want the scatter plot go in
    type (str): distance0.5, sasaLig or currentEnergy - different possibilities for the scatter plot
    """
    for key, value in data_dict.items():
        if "original" not in key:
            pele_profile_single(data_dict["original"], key, types, name, value, dpi)


def all_profiles(data_dict, name, dpi=1000):
    """
    Creates all the possible scatter plots for the same mutated position
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    name (str): name for the folders where you want the scatter plot go in
    """
    types = ["distance0.5", "sasaLig", "currentEnergy"]
    for x in types:
        pele_profiles(data_dict, name, x, dpi)


def extract_snapshot_from_pdb(simulation_folder, f_id, output, mutation, step, dist, out_freq=1):
    """
    Extracts PDB files from trajectories
    simulation_folder (str): Path to the simulation folder
    f_id (str): trajectory file ID
    output (str): The folder name for the results of the different simulations
    step (int): The step in the trajectory you want to keep
    out_freq (int): How frequent the steps are saved, in PELE every 1 step is saved
    mutation (str): The folder name for the results of one of the simulations
    dist (float): The distance between ligand and protein (used as name for the result file - not essential)
    """
    if not os.path.exists("results/distances_{}/{}_pdbs".format(output, mutation)):
        os.makedirs("results/distances_{}/{}_pdbs".format(output, mutation))

    f_in = glob("{}/output/0/*trajectory*_{}.*".format(simulation_folder, f_id))
    if len(f_in) == 0:
        sys.exit("Trajectory_{} not found. Be aware that PELE trajectories must contain the label 'trajectory' in "
                 "their file name to be detected".format(f_id))
    f_in = f_in[0]
    with open(f_in, 'r') as input_file:
        file_content = input_file.read()
    trajectory_selected = re.search(r'MODEL\s+{}(.*?)ENDMDL'.format(int((step/out_freq)+1)), file_content, re.DOTALL)

    # Output Snapshot
    traj = []
    path_ = "results/distances_{}/{}_pdbs".format(output, mutation)
    with open(os.path.join(path_, "traj{}_step{}_dist{}.pdb".format(f_id, step, round(dist, 2))), 'w') as f:
        traj.append("MODEL     {}".format(int((step/out_freq)+1)))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))


def extract_10_pdb_single(data, simulation_folder, output, mutation, out_freq=1):
    """
    Extracts the top 10 distances from one simulation
    data (SimulationData): A simulationData object that holds information of the simulation
    simulation_folder (str): Path to the simulation folders
    output (str): Folder name to store the results from different simulations
    mutation (str): Name for the folder to store results for one of the simulations
    out_freq (int): How frequent the steps are saved, in PELE every 1 step is saved
    """
    for ind in data.trajectory.index:
        ids = data.trajectory["ID"][ind]
        step = int(data.trajectory["numberOfAcceptedPeleSteps"][ind])
        dist = data.trajectory["distance0.5"][ind]
        extract_snapshot_from_pdb(simulation_folder, ids, output, mutation, step, dist, out_freq=out_freq)


def extract_10_pdbs_staturated(data_dict, folders, out_freq=1):
    """
    Extracts the top 10 distances for every 19 mutations at the same position
    data_dict (dict): A dictionary that contains SimulationData objects from the 19 simulation folders
    folders (str): Folder that has the results from different simulations at the same position
    out_freq (int): How frequent the steps are saved, in PELE every 1 step is saved
    """
    for folder in glob("{}/PELE_*".format(folders)):
        name = basename(folder)[5:]
        output = folder.split("/")[0]
        extract_10_pdb_single(data_dict[name], folder, output, mutation=name, out_freq=out_freq)


def consecutive_analysis(file_name, dpi=1000, out_freq=1):
    """
    Creates all the plots for the different mutated positions
    file_name (str): A file that contains the names of the different folders where the PELE simulation folders are in
    """
    if os.path.exists(file_name):
        with open("{}".format(file_name), "r") as pele:
            pele_folders = pele.readlines()
        for folders in pele_folders:
            folders = folders.strip("\n")
            data_dict = analyse_all(folders)
            box_plot(data_dict, folders, dpi)
            all_profiles(data_dict, folders, dpi)
            extract_10_pdbs_staturated(data_dict, folders, out_freq)
    else:
        raise OSError("No file {}".format(file_name))


def main():
    folder, dpi = parse_args()
    consecutive_analysis(folder, dpi)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
