from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--pele", required=False, default="./",
                        help="Include the folder path where the pele simulations are located")
    args = parser.parse_args()

    return args.pele


class SimulationData:
    def __init__(self, folder):
        """folder(str): path to the simulation folder"""
        self.folder = folder
        self.dataframe = None
        self.distance = None
        self.distribution = None

    def filtering(self):
        reports = []
        for files in glob("{}/output/0/report_*".format(self.folder)):
            reports.append(pd.read_csv(files, sep="    ", engine="python"))

        self.dataframe = pd.concat(reports)
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.head(len(self.dataframe)*20/100)
        self.distance = self.dataframe["distance0.5"].copy()
        self.distance.sort_values(inplace=True)
        self.distance.reset_index(drop=True, inplace=True)
        self.distance = self.distance.head(min(100, len(self.distance)))

        if "original" in self.folder:
            self.distance = self.distance[0].copy()

        return self.dataframe, self.distance

    def set_distribution(self, num):
        self.distribution = self.distance - num


def analyse_all(folders="."):
    """folders (str): path to the different PELE simulation folders to be analyzed"""
    data_dict = {}
    for folder in glob("{}/PELE_*".format(folders)):
        name = basename(folder)
        data = SimulationData(folder)
        data.filtering()
        data_dict[name[5:]] = data
        if "original" in folder:
            distance = data.distance
        else:
            distance = 3

    plot_dict = {}
    for key, value in data_dict.items():
        if "original" not in key:
            value.set_distribution(distance)
            plot_dict[key] = value.distribution

    plot_dataframe = pd.DataFrame(plot_dict)

    return data_dict, plot_dataframe


def box_plot(dataframe):
    """dataframe (tabular data): Any tabular structure that corresponds to the wide-form data accepted by seaborn"""
    sns.set(font_scale=1.8)
    sns.set_style("white")
    sns.set_context("paper")
    ax = sns.catplot(data=dataframe, kind="box", palette="Set2", height=4.5, aspect=2.3)
    ax.set(title="Distance variantion with respect to original")
    ax.set_ylabels("Distance", fontsize=8)
    ax.set_xlabels("Mutations", fontsize=8)
    ax.set_xticklabels(fontsize=6)
    ax.set_yticklabels(fontsize=6)
    ax.savefig("distance.png", dpi=200)

    return ax


def main():
    """folders (str): path to the different PELE simulation folders"""
    folder = parse_args()
    data_dict, plot_dataframe = analyse_all(folder)
    ax = box_plot(plot_dataframe)

    return data_dict, plot_dataframe, ax


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
