import argparse
from mutate_pdb import generate_mutations, generate_multiple_mutations
from pele_files import create_20sbatch
from subprocess import call
from os.path import abspath, basename
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, nargs="+",
                        help="Include one or more chain IDs and positions -> Chain ID:position")
    parser.add_argument("--ligchain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--ligname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=24, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--test", required=False, action="store_true")
    parser.add_argument("--cu", required=False, action="store_true")
    parser.add_argument("--multiple", required=False, action="store_true")
    parser.add_argument("--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")

    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atom1, args.atom2, args.cpus, args.test,
            args.cu, args.multiple, args.seed]


def submit(slurm_folder):
    """
    Given a folder submits the job to the supercomputer
    """
    for files in slurm_folder:
        call(["sbatch", "{}".format(files)])


def side_function(input_):
    """
    input_ (str): The wild type PDB file path
    Put all the necessary previous steps here
    """
    input_ = abspath(input_)
    base = basename(input_)
    base = base.replace(".pdb", "")
    if not os.path.exists("mutations_{}".format(base)):
        os.mkdir("mutations_{}".format(base))
    os.chdir("mutations_{}".format(base))

    return input_


def pele_folders(input_, file_list):
    """
    Creates a file with the names of the different folders where the pele simulations are contained
    input_ (str): The wild type PDB file pat
    file_list (list): list of pdb files created during the saturated mutagenesis
    """
    os.chdir("../")
    base = basename(input_)
    base = base.replace(".pdb", "")
    count = 0
    folder = []
    for files in file_list:
        name = basename(files)
        name = name.replace(".pdb", "")
        if not count:
            hold = "bla"
            count += 1
        if name != "original" and hold != name[:-1]:
            hold = name[:-1]
            folder.append("mutations_{}/{}\n".format(base, hold))
    with open("dirnames_{}.txt".format(base), "w") as txt:
        txt.writelines(folder)


def main():
    input_, position, ligchain, ligname, atom1, atom2, cpus, test, cu, multiple, seed = parse_args()
    input_ = side_function(input_)
    if multiple:
        pdb_names = generate_multiple_mutations(input_, position, hydrogens=True)
    else:
        pdb_names = generate_mutations(input_, position, hydrogens=True)

    slurm_files = create_20sbatch(ligchain, ligname, atom1, atom2, cpus=cpus, test=test, initial=input_,
                                  file_list=pdb_names, cu=cu, seed=seed)
    submit(slurm_files)
    pele_folders(input_, pdb_names)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
