"""
This script is designed to run saturated_mutagenesis through the command-line.
"""

__author__ = "Ruite Xiang"
__license__ = "MIT"
__maintainer__ = "Ruite Xiang"
__email__ = "ruite.xiang@bsc.es"


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
    parser.add_argument("--hydrogen", required=False, action="store_false", help="leave it to default")
    parser.add_argument("--consec", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")

    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atom1, args.atom2, args.cpus, args.test,
            args.cu, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen, args.consec]


class SimulationRunner:
    """
    A class that configures and runs simulations
    """

    def __init__(self, input_, dir_=None, single=None, nord=False):
        """
        Initialize the Simulation Runner class
        Parameters
        ___________
        input_: str
            The path to the PDB file
        dir_: str, optional
            The name of the directory for the simulations to run and the outputs to be stored
        nord: bool, optional
            Set to True if you want to run the simulation on NORDIII
        """

        self.input = input_
        self.dir = dir_
        self.single = single
        self.nord = nord

    def side_function(self):
        """
        Change the working directory to store all the simulations in one place
        Returns
        _______
        input_: str
            The new path of the input
        """
        self.input = abspath(self.input)
        if not self.dir:
            base = basename(self.input)
            base = base.replace(".pdb", "")
        else:
            base = self.dir
        if not os.path.exists("{}_mutations".format(base)):
            os.mkdir("{}_mutations".format(base))
        os.chdir("{}_mutations".format(base))

        return self.input

    def pele_folders(self, pdb_list):
        """
        Creates a file with the names of the different folders where the pele simulations are contained
        Parameters
        ___________
        pdb_list: list[path]
            list of pdb files path created during the saturated mutagenesis
        single: str
            Anything that indiucates that the plurizymes is used
        """
        os.chdir("../")
        if not self.dir:
            base = basename(self.input)
            base = base.replace(".pdb", "")
        else:
            base = basename(self.dir)
        hold = "bla"
        folder = []
        if not self.single:
            for files in pdb_list:
                name = basename(files).replace(".pdb", "")
                if name != "original" and hold != name[:-1]:
                    hold = name[:-1]
                    folder.append("{}_mutations/{}\n".format(base, hold))
            dirname = "dirnames_{}.txt".format(base)
            with open(dirname, "w") as txt:
                txt.writelines(folder)

            return dirname

    def submit(self, slurm_folder):
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
            if not self.nord:
                call(["sbatch", "{}".format(files)])
            else:
                os.system("bsub < {}".format(files))


def main():
    input_, position, ligchain, ligname, atom1, atom2, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, \
    hydrogen, consec = parse_args()
    simulation = SimulationRunner(input_, dir_, nord=nord)
    input_ = simulation.side_function()
    pdb_names = generate_mutations(input_, position, hydrogens=hydrogen, multiple=multiple, folder=pdb_dir, consec=consec)
    slurm_files = create_20sbatch(ligchain, ligname, atom1, atom2, cpus=cpus, test=test, initial=input_,
                                  file_=pdb_names, cu=cu, seed=seed, nord=nord)
    simulation.submit(slurm_files)
    simulation.pele_folders(pdb_names)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
