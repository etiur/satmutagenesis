import argparse
from mutate_pdb import generate_mutations
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
    parser.add_argument("--cu", required=False, action="store_true", help="used if there are copper in the system")
    parser.add_argument("--test", required=False, action="store_true", help="Used if you want to run a test before")
    parser.add_argument("--nord", required=False, action="store_true",
                        help="used if LSF is the utility managing the jobs")
    parser.add_argument("--multiple", required=False, action="store_true",
                        help="if you want to mutate 2 residue in the same pdb")
    parser.add_argument("--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    parser.add_argument("--dir", required=False,
                        help="The name of the folder for all the simulations")
    parser.add_argument("--pdb_dir", required=False, default="pdb_files",
                        help="The name for the mutated pdb folder")
    parser.add_argument("--hydrogen", required=False, action="store_false")
    parser.add_argument("--consec", required=False, action="store_true")

    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atom1, args.atom2, args.cpus, args.test,
            args.cu, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen, args.consec]


def submit(slurm_folder, nord=False):
    """
    Given a folder submits the job to the supercomputer

    Parameters
    __________
    slurm_folder: list[path]
        A list of the slurm files path's
    nord: bool, optional
        True if it will run on NORD
    """
    for files in slurm_folder:
        if not nord:
            call(["sbatch", "{}".format(files)])
        else:
            os.system("bsub < {}".format(files))


def side_function(input_, dir_=None):
    """
    Put all the necessary previous steps here

    Parameters
    ___________
    input_: str
        The wild type PDB file path
    dir_: str, optional
        Name of the folder ofr the simulations

    Returns
    _______
    input_: str
        The new path of the input
    """
    input_ = abspath(input_)
    if not dir_:
        base = basename(input_)
        base = base.replace(".pdb", "")
    else:
        base = dir_
    if not os.path.exists("mutations_{}".format(base)):
        os.mkdir("mutations_{}".format(base))
    os.chdir("mutations_{}".format(base))

    return input_


def pele_folders(input_, file_list, dir_=None):
    """
    Creates a file with the names of the different folders where the pele simulations are contained

    Parameters
    ___________
    input_: str
        The wild type PDB file path
    file_list: list[path]
        list of pdb files path created during the saturated mutagenesis
    dir_: str, optional
        Name of the folder ofr the simulations
    """
    os.chdir("../")
    if not dir_:
        base = basename(input_)
        base = base.replace(".pdb", "")
    else:
        base = dir_
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
    input_, position, ligchain, ligname, atom1, atom2, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, \
    hydrogen, consec = parse_args()
    input_ = side_function(input_, dir_)
    pdb_names = generate_mutations(input_, position, hydrogens=hydrogen, multiple=multiple, folder=pdb_dir, consec=consec)
    slurm_files = create_20sbatch(ligchain, ligname, atom1, atom2, cpus=cpus, test=test, initial=input_,
                                  file_=pdb_names, cu=cu, seed=seed, nord=nord)
    submit(slurm_files, nord)
    pele_folders(input_, pdb_names, dir_)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
