from glob import glob
import pandas as pd
import seaborn as sns
import argparse
from os.path import basename, dirname, abspath
import os
import matplotlib.pyplot as plt
import sys
import re
from fpdf import FPDF
import logging


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("--pele", required=True,
                        help="Include a file with names of the different folders with PELE simulations inside")
    parser.add_argument("--dpi", required=False, default=1000, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("--distance", required=False, default=30, type=int,
                        help="Set how many data points are used for the boxplot")
    parser.add_argument("--trajectory", required=False, default=10, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--out", required=False, default="summary",
                        help="Name of the summary file created at the end of the analysis")
    parser.add_argument("--folder", required=False,
                        help="Name of the results folder")
    args = parser.parse_args()

    return args.pele, args.dpi, args.distance, args.trajectory, args.out, args.folder


class SimulationData:
    def __init__(self, folder, points=30, pdb=10):
        """
        folder (str):  path to the simulation folder
        points (int): Number of points to consider for the boxplots
        pdb (int): how many pdbs to extract from the trajectories
        """
        self.folder = folder
        self.dataframe = None
        self.distance = None
        self.distribution = None
        self.profile = None
        self.trajectory = None
        self.points = points
        self.pdb = pdb
        self.binding = None
        self.bind_diff = None
        self.data = None

    def filtering(self):
        """
        Constructs a dataframe from all the reports in a PELE simulation folder with the best 20% binding energies
        and a Series with the 100 best ligand distances
        """
        pd.options.mode.chained_assignment = None
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
        self.trajectory = self.trajectory.iloc[:self.pdb]

        # For the box plots
        data_20 = self.dataframe.iloc[:len(self.dataframe) * 20 / 100]
        data_20.sort_values(by="distance0.5", inplace=True)
        data_20.reset_index(drop=True, inplace=True)
        data_20 = data_20.iloc[:min(self.points, len(data_20))]
        self.distance = data_20["distance0.5"].copy()
        self.binding = data_20["Binding Energy"].copy()
        self.binding.sort_values(inplace=True)
        self.binding.reset_index(drop=True, inplace=True)

        if "original" in self.folder:
            self.distance = self.distance.iloc[0]
            self.binding = self.binding.iloc[0]

    def set_distribution(self, original_distance):
        self.distribution = self.distance - original_distance

    def set_binding(self, original_binding):
        self.bind_diff = self.binding - original_binding


def analyse_all(folders=".", distance=30, trajectory=10):
    """
    folders (str): path to the different PELE simulation folders to be analyzed
    """
    data_dict = {}
    if len(folders.split("/")) > 1:
        mutation_dir = dirname(folders)
        original = SimulationData("{}/PELE_original".format(mutation_dir), points=distance, pdb=trajectory)
    else:
        original = SimulationData("PELE_original")
    original.filtering()
    data_dict["original"] = original
    for folder in glob("{}/PELE_*".format(folders)):
        name = basename(folder)
        data = SimulationData(folder, points=distance, pdb=trajectory)
        data.filtering()
        data.set_distribution(original.distance)
        data.set_binding(original.binding)
        data_dict[name[5:]] = data

    return data_dict


def box_plot(input_, data_dict, position_num, dpi=1000):
    """
    Creates a box plot of the 19 mutations from the same position
    input_ (str): name of the results folder
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    position_num (str): Position at the which the mutations occurred
    """
    if not os.path.exists("results_{}/Plots/box".format(input_)):
        os.makedirs("results_{}/Plots/box".format(input_))
    # create a dataframe with only the distance differences for each simulation
    plt.ioff()
    plot_dict_dist = {}
    plot_dict_bind = {}
    for key, value in data_dict.items():
        if "original" not in key:
            plot_dict_dist[key] = value.distribution
            plot_dict_bind[key] = value.bind_diff

    data_dist = pd.DataFrame(plot_dict_dist)
    data_bind = pd.DataFrame(plot_dict_bind)

    sns.set(font_scale=1.9)
    sns.set_style("ticks")
    sns.set_context("paper")
    # Distance boxplot
    ax = sns.catplot(data=data_dist, kind="box", palette="Accent", height=4.5, aspect=2.3)
    ax.set(title="{} distance variation with respect to wild type".format(position_num))
    ax.set_ylabels("Distance variation", fontsize=9)
    ax.set_xlabels("Mutations {}".format(position_num), fontsize=9)
    ax.set_xticklabels(fontsize=7)
    ax.set_yticklabels(fontsize=7)
    ax.savefig("results_{}/Plots/box/{}_distance.png".format(input_ ,position_num), dpi=dpi)
    # Binding energy Box plot
    ex = sns.catplot(data=data_bind, kind="box", palette="Accent", height=4.5, aspect=2.3)
    ex.set(title="{} Binding energy variation with respect to wild type".format(position_num))
    ex.set_ylabels("Binding energy variation", fontsize=9)
    ex.set_xlabels("Mutations {}".format(position_num), fontsize=9)
    ex.set_xticklabels(fontsize=7)
    ex.set_yticklabels(fontsize=7)
    ex.savefig("results_{}/Plots/box/{}_binding.png".format(input_, position_num), dpi=dpi)


def pele_profile_single(input_, wild, key, types, position_num, mutation, dpi=1000):
    """
    Creates a plot for a single mutation
    input_ (str): name of the results folder
    wild (SimulationData): SimulationData object that stores data for the wild type protein
    key (str): name of the mutation
    types (str): Type of scatter plot - distance0.5, sasaLig or currentEnergy
    position_num (str): name for the folder to keep the images from the different mutations
    mutation (SimulationData): SimulationData object that stores data for the mutated protein
    """
    # Configuring the plot
    plt.ioff()
    sns.set(font_scale=1.3)
    sns.set_style("ticks")
    sns.set_context("paper")
    original = wild.profile
    original.index = ["Wild type"] * len(original)
    distance = mutation.profile
    distance.index = [key] * len(distance)
    cat = pd.concat([original, distance], axis=0)
    cat.index.name = "Type"
    cat.reset_index(inplace=True)
    # Creating the scatter plots
    if types == "currentEnergy":
        if not os.path.exists("results_{}/Plots/scatter_{}_{}/{}".format(input_, position_num, types, "distance0.5")):
            os.makedirs("results_{}/Plots/scatter_{}_{}/{}".format(input_, position_num, types, "distance0.5"))
        if not os.path.exists("results_{}/Plots/scatter_{}_{}/{}".format(input_, position_num, types, "sasaLig")):
            os.makedirs("results_{}/Plots/scatter_{}_{}/{}".format(input_, position_num, types, "sasaLig"))

        norm = plt.Normalize(cat["sasaLig"].min(), cat["sasaLig"].max())
        norm2 = plt.Normalize(cat["distance0.5"].min(), cat["distance0.5"].max())
        ax = sns.relplot(x=types, y='Binding Energy', hue="sasaLig", style="Type", palette='RdBu', data=cat,
                         height=3.8, aspect=1.8, hue_norm=norm, s=100, linewidth=0)
        ex = sns.relplot(x=types, y='Binding Energy', hue="distance0.5", style="Type", palette='RdBu', data=cat,
                         height=3.8, aspect=1.8, hue_norm=norm2, s=100, linewidth=0)
        ex.set(title="{} scatter plot of binding energy vs {} ".format(key, types))
        ex.savefig("results_{}/Plots/scatter_{}_{}/{}/{}_{}.png".format(input_, position_num, types, "distance0.5", key, types), dpi=dpi)
        ax.savefig("results_{}/Plots/scatter_{}_{}/{}/{}_{}.png".format(input_, position_num, types, "sasaLig", key, types), dpi=dpi)

    else:
        if not os.path.exists("results_{}/Plots/scatter_{}_{}".format(input_, position_num, types)):
            os.makedirs("results_{}/Plots/scatter_{}_{}".format(input_, position_num, types))
        ax = sns.relplot(x=types, y='Binding Energy', hue="Type", style="Type", palette="Set1", data=cat,
                         height=3.8, aspect=1.8, s=100, linewidth=0)
        ax.set(title="{} scatter plot of binding energy vs {} ".format(key, types))
        ax.savefig("results_{}/Plots/scatter_{}_{}/{}_{}.png".format(input_, position_num, types, key, types), dpi=dpi)


def pele_profiles(input_, data_dict, position_num, types, dpi=1000):
    """
    Creates a scatter plot for each of the 19 mutations from the same position by comparing it to the wild type
    input_ (str): Name of the results folder
    data_dict (dict): A dictionary that contains SimulationData objects from the 19 simulation folders
    position_num (str): name for the folders where you want the scatter plot go in
    type (str): distance0.5, sasaLig or currentEnergy - different possibilities for the scatter plot
    """
    for key, value in data_dict.items():
        if "original" not in key:
            pele_profile_single(input_, data_dict["original"], key, types, position_num, value, dpi)


def all_profiles(input_, data_dict, position_num, dpi=1000):
    """
    Creates all the possible scatter plots for the same mutated position
    input_ (str): Name of the results folder
    data_dict (dict): A dictionary that contains SimulationData objects from the simulation folders
    position_num (str): name for the folders where you want the scatter plot go in
    """
    types = ["distance0.5", "sasaLig", "currentEnergy"]
    for x in types:
        pele_profiles(input_, data_dict, position_num, x, dpi)


def extract_snapshot_from_pdb(input_, simulation_folder, f_id, position_num, mutation, step, dist, bind):
    """
    Extracts PDB files from trajectories
    simulation_folder (str): Path to the simulation folder
    f_id (str): trajectory file ID
    position_num (str): The folder name for the output of this function for the different simulations
    step (int): The step in the trajectory you want to keep
    mutation (str): The folder name for the output of this function for one of the simulations
    dist (float): The distance between ligand and protein (used as name for the result file - not essential)
    bind (float): The binding energy between ligand and protein (used as name for the result file - not essential)
    """
    if not os.path.exists("results_{}/distances_{}/{}_pdbs".format(input_, position_num, mutation)):
        os.makedirs("results_{}/distances_{}/{}_pdbs".format(input_, position_num, mutation))

    f_in = glob("{}/output/0/*trajectory*_{}.*".format(simulation_folder, f_id))
    if len(f_in) == 0:
        sys.exit("Trajectory_{} not found. Be aware that PELE trajectories must contain the label 'trajectory' in "
                 "their file name to be detected".format(f_id))
    f_in = f_in[0]
    with open(f_in, 'r') as input_file:
        file_content = input_file.read()
    trajectory_selected = re.search(r'MODEL\s+{}(.*?)ENDMDL'.format(int(step)+1), file_content, re.DOTALL)

    # Output Snapshot
    traj = []
    path_ = "results_{}/distances_{}/{}_pdbs".format(input_, position_num, mutation)
    name = "traj{}_step{}_dist{}_bind{}.pdb".format(f_id, step, round(dist, 2), round(bind, 2))
    with open(os.path.join(path_, name), 'w') as f:
        traj.append("MODEL     {}".format(int(step)+1))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))


def extract_10_pdb_single(input_, data, simulation_folder, position_num, mutation):
    """
    Extracts the top 10 distances for one mutation
    input_ (str): Name of the results folder
    data (SimulationData): A simulationData object that holds information of the simulation
    simulation_folder (str): Path to the simulation folders
    position_num (str): The folder name for the output of this function for the different simulations
    mutation (str): Name for the folder to store results for one of the simulations
    """
    for ind in data.trajectory.index:
        ids = data.trajectory["ID"][ind]
        step = data.trajectory["numberOfAcceptedPeleSteps"][ind]
        dist = data.trajectory["distance0.5"][ind]
        bind = data.trajectory["Binding Energy"][ind]
        extract_snapshot_from_pdb(input_, simulation_folder, ids, position_num, mutation, step, dist, bind)


def extract_all(input_, data_dict, position_num):
    """
    Extracts the top 10 distances for the 19 mutations at the same position
    input_ (str): name of the results folder
    data_dict (dict): A dictionary that contains SimulationData objects from the 19 simulation folders
    position_num (str): Folder that has the results from the different simulations from the same residue number
    """
    for folder in glob("{}/PELE_*".format(position_num)):
        name = basename(folder)[5:]
        if len(folder.split("/")) > 2:
            output = folder.split("/")[1]
        else:
            output = folder.split("/")[0]
        extract_10_pdb_single(input_, data_dict[name], folder, output, mutation=name)


def create_report(input_, mutation, position_num, output="summary"):
    """
    Create pdf files with the plots of chosen mutations and the path to the
    input_ (str): Name of the results folder
    mutation (dict): {mutations: [distances, binding energies]}
    position_num (str): part of the path to the plots
    output (str): The pdf filename without the extension
    return: pdf file
    """
    pdf = FPDF()
    pdf.set_top_margin(17.0)
    pdf.set_left_margin(15.0)
    pdf.set_right_margin(15.0)
    pdf.add_page()
    # Title
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, "Best mutations in terms of distance and binding energy", align='C', ln=1)
    pdf.set_font('Arial', '', size=10)
    for mut, key in mutation.items():
        dis = key.distribution.median()
        bind = key.bind_diff.median()
        message = 'Mutation {}: median distance increment {}, median binding energy increment {}'.format(mut, dis, bind)
        pdf.ln(3)  # linebreaks
        pdf.cell(0, 5, message, ln=1)
    pdf.ln(8)  # linebreaks

    # Plots
    pdf.set_font('Arial', 'B', size=12)
    pdf.cell(0, 10, "Plots", align='C', ln=1)
    pdf.set_font('Arial', '', size=10)
    for mut, key in mutation.items():
        pdf.ln(3)
        pdf.cell(0, 10, "Plots {}".format(mut), ln=1)
        pdf.ln(3)
        plot1 = "results_{}/Plots/scatter_{}_{}/{}_{}.png".format(input_, position_num, "distance0.5", mut, "distance0.5")
        plot2 = "results_{}/Plots/scatter_{}_{}/{}_{}.png".format(input_, position_num, "sasaLig", mut, "sasaLig")
        plot3 = "results_{}/Plots/scatter_{}_{}/{}/{}_{}.png".format(input_, position_num, "currentEnergy", "sasaLig", mut, "currentEnergy")
        plot4 = "results_{}/Plots/scatter_{}_{}/{}/{}_{}.png".format(input_, position_num, "currentEnergy", "distance0.5", mut, "currentEnergy")
        pdf.image(plot1, w=180)
        pdf.ln(3)
        pdf.image(plot2, w=180)
        pdf.ln(1000000)  # page break
        pdf.ln(3)
        pdf.image(plot3, w=180)
        pdf.ln(3)
        pdf.image(plot4, w=180)
        pdf.ln(1000000)  # page break

    # Top poses
    for mut, key in mutation.items():
        pdf.set_font('Arial', 'B', size=12)
        pdf.cell(0, 10, "Path for the top poses", align='C', ln=1)
        pdf.set_font('Arial', size=10)
        pdf.ln(5)
        path = "results_{}/distances_{}/{}_pdbs".format(input_, position_num, mut)
        pdf.cell(0, 10, "Top poses {}: {} ".format(mut, abspath(path)), ln=1)
        pdf.ln(5)

    # Output report
    pdf.output("results_{}/{}_{}.pdf".format(input_, output, position_num), 'F')
    return output


def find_top_mutations(input_, data_dict, position_num, output="summary"):
    """
    Finds those mutations that decreases the binding distance and binding energy and create a report
    input_ (str): Name of the results folder
    data_dict (dict): A dictionary of SimulationData objects that holds information for all mutations
    position_num (str): Part of the path to the plots included at the reports
    output (str): Name of the reports created
    """
    # Find top mutations
    count = 0
    mutation_dict = {}
    for key, value in data_dict.items():
        if value.distribution.median() < 0 and value.bind_diff.median() < 0:
            mutation_dict[key] = value
            count += 1
    # Create a summary report with the top mutations
    if len(mutation_dict) != 0:
        logging.info("{} mutations at position {} estimated to improve the system".format(count, position_num))
        create_report(input_, mutation_dict, position_num, output)
    else:
        logging.warning("No residues at position {} estimated to improve the system".format(position_num))


def consecutive_analysis(file_name, dpi=1000, distance=30, trajectory=10, output="summary", input_=None):
    """
    Creates all the plots for the different mutated positions
    input_ (str): Name for the results folder
    file_name (str): A file that contains the names of the different folders where the PELE simulation folders are in
    dpi (int): The quality of the plots
    distance (int): how many points are used for the boxplots
    trajectory (int): how many top pdbs are extracted from the trajectories
    output (str): name of the output file for the pdfs
    """
    if os.path.exists(file_name):
        with open("{}".format(file_name), "r") as pele:
            pele_folders = pele.readlines()
        if not input_:
            input_ = pele_folders[0].strip("\n")
            input_ = basename(dirname(input_)).replace("mutations_", "")
        logging.basicConfig(filename='results_{}/top_mutations.log'.format(input_), level=logging.DEBUG)
        for folders in pele_folders:
            folders = folders.strip("\n")
            data_dict = analyse_all(folders, distance=distance, trajectory=trajectory)
            box_plot(input_, data_dict, folders, dpi)
            all_profiles(input_, data_dict, folders, dpi)
            extract_all(input_, data_dict, folders)
            find_top_mutations(input_, data_dict, folders, output)
    else:
        raise OSError("No file named {}".format(file_name))


def main():
    pele, dpi, distance, trajectory, out, folder = parse_args()
    consecutive_analysis(pele, dpi, distance, trajectory, out, folder)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
