from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename
import os
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--pele", required=False, default="./",
                        help="Include a file with names of the different folders with PELE simulations inside")
    args = parser.parse_args()

    return args.pele


class SimulationData:
    original_distance = None

    def __init__(self, folder):
        """
        folder (str):  path to the simulation folder
        """
        self.folder = folder
        self.dataframe = None
        self.distance = None
        self.distribution = None
        self.profile = None

    def filtering(self):
        """
        Constructs a dataframe from all the reports in a PELE simulation folder with the best 20% binding energies
        and a Series with the 100 best ligand distances
        """
        reports = []
        for files in glob("{}/output/0/report_*".format(self.folder)):
            reports.append(pd.read_csv(files, sep="    ", engine="python"))
        self.dataframe = pd.concat(reports)
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)

        # for the PELE profiles
        self.profile = self.dataframe.drop(["Step", '#Task', 'numberOfAcceptedPeleSteps'], axis=1, inplace=True)
        self.profile = self.profile.iloc[:len(self.profile)-49]
        
        # For the box plots
        data_20 = self.dataframe.head(len(self.dataframe) * 20 / 100)
        dist20 = data_20["distance0.5"].copy()
        dist20.sort_values(inplace=True)
        dist20.reset_index(drop=True, inplace=True)
        self.distance = dist20.head(min(100, len(dist20)))
        if "original" in self.folder:
            self.distance = dist20[0].copy()

    def set_distribution(self):
        """
        Stores the difference between the mutated ligand distances and the original ligand distance
        """
        self.distribution = self.distance - self.original_distance


def analyse_all(folders="."):
    """
    folders (str): path to the different PELE simulation folders to be analyzed
    """
    # create a SimulationData object for each pele simulation
    data_dict = {}
    original = SimulationData("PELE_original")
    original.filtering()
    SimulationData.original_distance = original.distance
    data_dict["original"] = original
    for folder in glob("{}/PELE_*".format(folders)):
        name = basename(folder)
        data = SimulationData(folder)
        data.filtering()
        data.set_distribution()
        data_dict[name[5:]] = data

    return data_dict


def box_plot(data_dict, name):
    """
    Creates a box plot of the 19 mutations from the same position
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    """
    if not os.path.exists("Plots"):
        os.mkdir("Plots")
    if not os.path.exists("Plots/box"):
        os.mkdir("Plots/box")
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
    ax.savefig("Plots/box/{}.png".format(name), dpi=1500)


def pele_profiles(data_dict, name, types):
    """
    Creates a scatter plot for each of the 19 mutations from the same position by comparing it to the wild type
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    name (str): name for the folders where you want the scatter plot go in
    type (str): distance0.5, sasaLig or currentEnergy - different possibilities for the scatter plot
    """
    if not os.path.exists("Plots"):
        os.mkdir("Plots")
    if not os.path.exists("Plots/scatter_{}".format(name)):
        os.mkdir("Plots/scatter_{}".format(name))

    sns.set(font_scale=1)
    sns.set_style("ticks")
    sns.set_context("paper")
    original = data_dict["original"].profile
    original.index = ["Wild type"] * len(original)
    for key, value in data_dict.items():
        if "original" not in key:
            distance = value.profile
            distance.index = [key] * len(distance)
            cat = pd.concat([original, distance], axis=0)
            cat.index.name = "Type"
            cat.reset_index(inplace=True)
            if types == "currentEnergy":
                norm = plt.Normalize(cat["sasaLig"].min(), cat["sasaLig"].max())
                norm2 = plt.Normalize(cat["distance0.5"].min(), cat["distance0.5"].max())

                ax = sns.relplot(x=types, y='Binding Energy', hue="sasaLig", style="Type", palette='RdBu', data=cat,
                                 height=3.8, aspect=1.8, hue_norm=norm, s=120)
                ex = sns.relplot(x=types, y='Binding Energy', hue="distance0.5", style="Type", palette='RdBu', data=cat,
                                 height=3.8, aspect=1.8, hue_norm=norm2, s=120)
                ex.set(title="{} scatter plot of binding energy vs {} ".format(key, types))
                ex.savefig("Plots/scatter_{}/{}_{}_{}.png".format(name, key, types, "distance0.5"), dpi=1500)

            else:
                ax = sns.relplot(x=types, y='Binding Energy', hue="Type", style="Type", palette="Set1", data=cat,
                                 height=3.8, aspect=1.8, s=120)
            ax.set(title="{} scatter plot of binding energy vs {} ".format(key, types))
            ax.savefig("Plots/scatter_{}/{}_{}.png".format(name, key, types), dpi=1500)


def all_profiles(data_dict, name):
    """
    Creates all the possible scatter plots
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    name (str): name for the folders where you want the scatter plot go in
    """
    types = ["distance0.5", "sasaLig", "currentEnergy"]
    for x in types:
        pele_profiles(data_dict, name, x)


def consecutive_analysis(file_name):
    """
    file_name (str): A file that contains the names of the different folders where the PELE simulation folders are in
    """
    with open("{}".format(file_name), "r") as pele:
        pele_folders = pele.readlines()
    for folders in pele_folders:
        folders = folders.strip("\n")
        data_dict = analyse_all(folders)
        box_plot(data_dict, folders)
        all_profiles(data_dict, folders)


def main():
    folder = parse_args()
    consecutive_analysis(folder)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
