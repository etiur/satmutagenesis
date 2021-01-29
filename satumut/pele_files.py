"""
This script is used to generate the yaml files for pele platform
"""

import argparse
import os
from helper import map_atom_string, isiterable
from os.path import basename, join, isfile, isdir


def parse_args():
    parser = argparse.ArgumentParser(description="Generate running files for PELE")
    # main required arguments
    parser.add_argument("--folder", required=True,
                        help="An iterable of the path to different pdb files, a name of the folder or a file with the "
                             "path to the different pdb files")
    parser.add_argument("-lc","--ligchain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("-ln","--ligname", required=True, help="The ligand residue name")
    parser.add_argument("-at","--atoms", required=True, nargs="+",
                        help="Series of atoms of the residues to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cu", required=False, action="store_true", help="used if there are copper in the system")
    parser.add_argument("-t","--test", required=False, action="store_true", help="Used if you want to run a test before")
    parser.add_argument("-n","--nord", required=False, action="store_true",
                        help="used if LSF is the utility managing the jobs")
    parser.add_argument("-s","--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    parser.add_argument("-st", "--steps", required=False, type=int, default=700,
                        help="The number of PELE steps")
    args = parser.parse_args()

    return [args.folder, args.ligchain, args.ligname, args.atoms, args.cpus, args.test, args.cu,
            args.seed, args.nord, args.steps]


class CreateYamlFiles:
    """
    Creates the 2 necessary files for the pele simulations
    """

    def __init__(self, input_, ligchain, ligname, atoms, cpus=25,
                 test=False, initial=None, cu=False, seed=12345, nord=False, steps=700):
        """
        Initialize the CreateLaunchFiles object

        Parameters
        ___________
        input_: str
            A PDB file's path
        ligchain: str
            the chain ID where the ligand is located
        ligname: str
            the residue name of the ligand in the PDB
        atoms: list[str]
            list of atom of the residue to follow, in this format --> chain ID:position:atom name
        cpus: int, optional
            How many cpus do you want to use
        test: bool, optional
            Setting the simulation to test mode
        initial: file, optional
            The initial PDB file before the modification by pmx
        cu: bool, optional
            Set it to true if there are coppers in the system
        seed: int, optional
            A seed number to make the simulations reproducible
        nord: bool, optional
            True if the system is managed by LSF
        steps: int, optional
            The number of PELE steps
        """

        self.input = input_
        self.ligchain = ligchain
        self.ligname = ligname
        self.atoms = atoms
        self.cpus = cpus
        self.test = test
        self.yaml = None
        self.initial = initial
        self.cu = cu
        self.seed = seed
        self.nord = nord
        self.steps = steps

    def _match_dist(self):
        """
        match the user coordinates to pmx PDB coordinates
        """
        if self.initial:
            for i in range(len(self.atoms)):
                self.atoms[i] = map_atom_string(self.atoms[i], self.initial, self.input)
        else:
            pass

    def input_creation(self, name):
        """
        create the .yaml input files for PELE

        Parameters
        ___________
        yaml_name: str
            Name for the input file for the simulation
        """
        self._match_dist()

        if not os.path.exists("yaml_files"):
            os.mkdir("yaml_files")
        self.yaml = "yaml_files/{}.yaml".format(name)
        with open(self.yaml, "w") as inp:
            lines = ["system: '{}'\n".format(self.input), "chain: '{}'\n".format(self.ligchain),
                     "resname: '{}'\n".format(self.ligname), "induced_fit_exhaustive: true\n",
                     "seed: {}\n".format(self.seed), "steps: {}\n".format(self.steps), "atom_dist:\n"]
            lines_atoms = ["- '{}'\n".format(atoms) for atoms in self.atoms]
            lines.extend(lines_atoms)
            if not self.nord:
                lines.append("usesrun: true\n")
            if name != "original":
                lines.append("working_folder: {}/PELE_{}\n".format(name[:-1], name))
            else:
                lines.append("working_folder: PELE_{}\n".format(name))
            if self.test:
                lines.append("test: true\n")
                self.cpus = 5
            lines2 = ["cpus: {}\n".format(self.cpus),
                      "pele_license: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/license'\n",
                      "pele_exec: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/bin/PELE-1.6.1_mpi'\n"]
            if self.cu:
                path = "/gpfs/projects/bsc72/ruite/examples/cuz"
                lines2.append("templates:\n- '{}'\n".format(path))
            lines.extend(lines2)
            inp.writelines(lines)

        return self.yaml


def create_20sbatch(ligchain, ligname, atoms, file_, cpus=25, test=False, initial=None,
                    cu=False, seed=12345, nord=False, steps=700):
    """
    creates for each of the mutants the yaml and slurm files

    Parameters
    ___________
    ligchain: str
        the chain ID where the ligand is located
    ligname: str
        the residue name of the ligand in the PDB
    atoms: list[str]
        list of atom of the residue to follow, in this format --> chain ID:position:atom name
    file_: iterable (not string or dict), dir or a file
        An iterable of the path to different pdb files, a name of the folder
        or a file of the path to the different pdb files
    cpus: int, optional
        how many cpus do you want to use
    test: bool, optional
        Setting the simulation to test mode
    initial: file, optional
        The initial PDB file before the modification by pmx if the residue number are changed
    cu: bool, optional
        Set it to true if there are coppers in the system
    seed: int, optional
        A seed number to make the simulations reproducible
    nord: bool, optional
        True if the system is managed by LSF
    steps: int, optional
            The number of PELE steps

    Returns
    _______
    slurm_files: list[path]
        A list of the files generated
    """
    slurm_files = []
    if isdir(str(file_)):
        file_list = list(filter(lambda x: ".pdb" in x, os.listdir(file_)))
        file_list = [join(file_, files) for files in file_list]
    elif isfile(str(file_)):
        with open("{}".format(file_), "r") as pdb:
            file_list = pdb.readlines()
    elif isiterable(file_):
        file_list = file_[:]
    else:
        raise Exception("No directory or iterable passed")

    # Create the launching files
    yaml_files = []
    for files in file_list:
        files = files.strip("\n")
        name = basename(files).replace(".pdb", "")
        run = CreateYamlFiles(files, ligchain, ligname, atoms, cpus, test=test,
                              initial=initial, cu=cu, seed=seed, nord=nord, steps=steps)
        yaml = run.input_creation(name)
        yaml_files.append(yaml)

    return yaml_files


def main():
    folder, ligchain, ligname, atoms, cpus, test, cu, seed, nord, steps = parse_args()
    yaml_files = create_20sbatch(ligchain, ligname, atoms,
                                 cpus=cpus, file_=folder, test=test, cu=cu, seed=seed, nord=nord, steps=steps)

    return yaml_files


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    yaml_list, slurm_list = main()
