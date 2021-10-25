"""
This script is used to analyse the results of the simulations
"""

import glob
import argparse
from .PELEAddMetric import MetricCalculator


def parse_args():
    parser = argparse.ArgumentParser(description="Analyse the different PELE simulations and create plots")
    # main required arguments
    parser.add_argument("input_dir", help="Path to the folders with PELE simulations inside")
    parser.add_argument("round_kind",
                        help="Aminoacid to mutate, options are (SER, HIS, ASP, GLU)")
    parser.add_argument("plurizyme_at_and_res",
                        help="Specify the chain ID, residue number and the PDB atom name that"
                             "will set the list of the neighbouring residues for the"
                             "next round. Example: chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")

    args = parser.parse_args()

    return [args.input_dir, args.round_kind, args.plurizyme_at_and_res,
            args.cpus, args.xtc]

def PELE_PLURIZYMER_add_metrics(_dir, Round, plurizyme_at_and_res, xtc, total_cpus):

    if Round == "S" or Round == "SER":
        List_of_mutants = glob.glob("{}/round_1/Subset_*/output/*S".format(_dir))
        List_of_mutants.sort()

        for i, mutant in enumerate(List_of_mutants):
            metrics_round_1 = MetricCalculator(mutant,list_of_residue_numbers= [1, mutant.split("/")[-1][1:-1]], \
            list_of_atom_names=[plurizyme_at_and_res.split(":")[-1],"OG"],list_of_residue_names=[plurizyme_at_and_res.split(":")[1],"SER"], \
            column_name= "Serine-Substrate Distance", report_format="report*",trajectory_format= "xtc" if xtc else "pdb", topology=("/".join(mutant.split("/")[:-1])+"//topologies/topology_{}.pdb".format(i)), \
            output_name="metric_round_{}".format(Round), n_processors= total_cpus)
            metrics_round_1.AddLigMetric()


    elif Round == "H" or Round == "HIS":

        List_of_mutants = glob.glob("{}/round_2/*/Subset_*/output/*H".format(_dir))
        List_of_mutants.sort()

        for i, mutant in enumerate(List_of_mutants):
            # TODO: Get the the residue number of the previous round for the list_of_residue_numbers attribute of the MetricCalculator class
            metrics_round_2 = MetricCalculator(mutant,list_of_residue_numbers= [mutant.split("/")[-1].split("_")[-2][1:-1], mutant.split("/")[-1].split("_")[-1][1:-1]], \
            list_of_atom_names=["HG","NE2"],list_of_residue_names=["SER", "HIS"], \
            column_name= "Serine-Histidine Distance", report_format="report*",trajectory_format= "xtc" if xtc else "pdb", topology=("/".join(mutant.split("/")[:-1])+ "/topologies/topology_{}.pdb".format(i)), \
            output_name="metric_round_{}".format(Round), n_processors= total_cpus)
            metrics_round_2.AddMetric()

    elif Round == "D" or Round == "ASP":
        List_of_mutants = glob.glob("{}/round_3/*/Subset_*/output/*D".format(_dir))
        List_of_mutants.sort()

        for i, mutant in enumerate(List_of_mutants):
            # TODO: Get the the residue number of the previous round for the list_of_residue_numbers attribute of the MetricCalculator class
            print(mutant,[mutant.split("/")[-1].split("_")[-2][1:-1], mutant.split("/")[-1].split("_")[-1][1:-1]])
            metrics_round_31 = MetricCalculator(mutant,list_of_residue_numbers= [mutant.split("/")[-1].split("_")[-2][1:-1], mutant.split("/")[-1].split("_")[-1][1:-1]], \
            list_of_atom_names=["HD1","OD1"],list_of_residue_names=["HIS", "ASP"], \
            column_name= "Aspartate-Histidine Distance 1", report_format="report*",trajectory_format= "xtc" if xtc else "pdb", topology=("/".join(mutant.split("/")[:-1])+ "/topologies/topology_{}.pdb".format(i)), \
            output_name="metric_round_{}".format(Round), n_processors= total_cpus)
            metrics_round_31.AddMetric()
            metrics_round_32 = MetricCalculator(mutant,list_of_residue_numbers= [mutant.split("/")[-1].split("_")[-2][1:-1], mutant.split("/")[-1].split("_")[-1][1:-1]], \
            list_of_atom_names=["HD1","OD2"],list_of_residue_names=["HIS", "ASP"], \
            column_name= "Aspartate-Histidine Distance 2", report_format="report*metric_round*",trajectory_format= "xtc" if xtc else "pdb", topology=("/".join(mutant.split("/")[:-1])+ "/topologies/topology_{}.pdb".format(i)), \
            output_name="2", n_processors= total_cpus)
            metrics_round_32.AddMetric()

    elif Round == "E" or Round == "GLU":
        # TODO: The Subset_1 should be replaced by the list of mutants from the previous round
        List_of_mutants = glob.glob("{}/round_3/*/Subset_*/output/*E".format(_dir))
        List_of_mutants.sort()

        for i, mutant in enumerate(List_of_mutants):
            # TODO: Get the the residue number of the previous round for the list_of_residue_numbers attribute of the MetricCalculator class
            metrics_round_31 = MetricCalculator(mutant,list_of_residue_numbers= [mutant.split("/")[-1].split("_")[-1][1:-1], mutant.split("/")[-1][1:-1]], \
            list_of_atom_names=["HD1","OE1"],list_of_residue_names=["HIS", "GLU"], \
            column_name= "Glutamate-Histidine Distance 1", report_format="report*",trajectory_format= "xtc" if xtc else "pdb", topology=("/".join(mutant.split("/")[:-1])+ "/topologies/topology_{}.pdb".format(i)), \
            output_name="metric_round_{}".format(Round), n_processors= total_cpus)
            metrics_round_31.AddMetric()
            metrics_round_32 = MetricCalculator(mutant,list_of_residue_numbers= [mutant.split("/")[-1].split("_")[-1][1:-1], mutant.split("/")[-1][1:-1]], \
            list_of_atom_names=["HD1","OE2"],list_of_residue_names=["HIS", "GLU"], \
            column_name= "Glutamate-Histidine Distance 2", report_format="report*metric_round*",trajectory_format= "xtc" if xtc else "pdb", topology=("/".join(mutant.split("/")[:-1])+ "//topologies/topology_{}.pdb".format(i)), \
            output_name="2", n_processors= total_cpus)
            metrics_round_32.AddMetric()



def main():

    directory, round_kind, plurizyme_at_and_res, cpus, xtc = parse_args()
    PELE_PLURIZYMER_add_metrics(directory, round_kind, plurizyme_at_and_res, xtc, cpus)

if __name__ == "__main__":
    main()
