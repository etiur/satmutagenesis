from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename


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
        self.dataframe = self.dataframe.head(len(self.dataframe) * 20/100)
        self.distance = self.dataframe["distance0.5"].copy()
        self.distance.sort_values(inplace=True)
        self.distance.reset_index(drop=True, inplace=True)
        self.distance = self.distance.head(min(100, len(self.distance)))

        if "original" in self.folder:
            self.distance = self.distance[0].copy()

        return self.dataframe, self.distance

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
    for folder in glob("{}/PELE_*".format(folders)):
        name = basename(folder)
        data = SimulationData(folder)
        data.filtering()
        data.set_distribution()
        data_dict[name[5:]] = data

    # create a dataframe with only the distance differences
    plot_dict = {}
    for key, value in data_dict.items():
        if "original" not in key:
            plot_dict[key] = value.distribution
    plot_dataframe = pd.DataFrame(plot_dict)

    return data_dict, plot_dataframe


def box_plot(dataframe, name):
    """
    dataframe (tabular data): Any tabular structure that corresponds to the wide-form data accepted by seaborn
    """
    sns.set(font_scale=1.8)
    sns.set_style("white")
    sns.set_context("paper")
    ax = sns.catplot(data=dataframe, kind="box", palette="Set2", height=4.5, aspect=2.3)
    ax.set(title="Distance variantion with respect to original")
    ax.set_ylabels("Distance", fontsize=9)
    ax.set_xlabels("Mutations", fontsize=9)
    ax.set_xticklabels(fontsize=7)
    ax.set_yticklabels(fontsize=7)
    ax.savefig("{}.png".format(name), dpi=1500)


def consecutive_analysis(file_name):
    """
    file_name (str): A file that contains the names of the different folders where the PELE simulation folders are in
    """
    with open("{}".format(file_name), "r") as pele:
        pele_folders = pele.readlines()
    for folders in pele_folders:
        folders = folders.strip("\n")
        data_dict, plot_dataframe = analyse_all(folders)
        box_plot(plot_dataframe, folders)


def main():
    folder = parse_args()
    consecutive_analysis(folder)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
