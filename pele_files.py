import argparse
import os
from glob import glob
from helper import map_atom_string
from os.path import basename


def parse_args():
    parser = argparse.ArgumentParser(description="Generate running files for PELE")
    # main required arguments
    parser.add_argument("--folder", required=False, default="pdb_files",
                        help="Include the folder where the pdb files are located")
    parser.add_argument("--ligchain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--ligname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=24, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--cu", required=False, action="store_true")
    parser.add_argument("--test", required=False, action="store_true")
    args = parser.parse_args()

    return args.folder, args.ligchain, args.ligname, args.atom1, args.atom2, args.cpus, args.test, args.cu


class CreateLaunchFiles:
    def __init__(self, input_, ligchain, ligname, atom1, atom2, cpus=24, test=False, initial=None, cu=False):
        """
        input_ (str): PDB files path
        ligchain (str): the chain ID where the ligand is located
        ligname (str): the residue name of the ligand in the PDB
        atom1 (str): atom of the residue to follow in this format --> chain ID:position:atom name
        atom2 (str): atom of the ligand to follow in this format --> chain ID:position:atom name
        cpus (str or int): how many cpus do you want to use
        """
        self.input = input_
        self.ligchain = ligchain
        self.ligname = ligname
        self.atom1 = atom1
        self.atom2 = atom2
        self.cpus = cpus
        self.test = test
        self.yaml = None
        self.slurm = None
        self.initial = initial
        self.cu = cu

    def match_dist(self):
        """
        match the user coordinates to pmx PDB coordinates
        """
        if self.initial:
            self.atom1 = map_atom_string(self.atom1, self.initial, self.input)
            self.atom2 = map_atom_string(self.atom2, self.initial, self.input)
        else:
            pass

    def input_creation(self, yaml_name):
        """
        create the .yaml input files for PELE
        """
        if not os.path.exists("yaml_files"):
            os.mkdir("yaml_files")
        self.yaml = "yaml_files/{}.yaml".format(yaml_name)
        with open(self.yaml, "w") as inp:
            lines = ["system: '{}'\n".format(self.input), "chain: '{}'\n".format(self.ligchain),
                     "resname: '{}'\n".format(self.ligname), "induced_fit_exhaustive: true\n", "seed: 12345\n",
                     "usesrun: true\n"]
            inp.writelines(lines)
            if yaml_name != "original":
                inp.write("working_folder: {}/PELE_{}\n".format(yaml_name[:-1], yaml_name))
            else:
                inp.write("working_folder: PELE_{}\n".format(yaml_name))
            if self.test:
                inp.write("test: true\n")
                self.cpus = 5
            lines2 = ["cpus: {}\n".format(self.cpus), "atom_dist:\n- '{}'\n- '{}'\n".format(self.atom1, self.atom2),
                      "pele_license: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/license'\n",
                      "pele_exec: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/bin/PELE-1.6.1_mpi'\n"]
            if self.cu:
                path = "/gpfs/projects/bsc72/ruite/examples/cuz"
                inp.write("templates:\n- '{}'\n".format(path))
            inp.writelines(lines2)

    def slurm_creation(self, slurm_name):
        """
        Creates the slurm running files for PELE
        """
        if not os.path.exists("slurm_files"):
            os.mkdir("slurm_files")
        self.slurm = "slurm_files/{}.sh".format(slurm_name)
        with open(self.slurm, "w") as slurm:
            lines = ["#!/bin/bash\n", "#SBATCH -J PELE\n", "#SBATCH --output={}.out\n".format(slurm_name),
                     "#SBATCH --error={}.err\n".format(slurm_name)]
            slurm.writelines(lines)
            if self.test:
                slurm.write("#SBATCH --qos=debug\n")
                self.cpus = 5
            lines2 = ["#SBATCH --ntasks={}\n\n".format(self.cpus), 'module purge\n',
                      'export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"\n',
                      'export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"\n',
                      'export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH\n', 'module load impi\n',
                      'module load intel mkl impi gcc # 2> /dev/null\n', 'module load boost/1.64.0\n',
                      '/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8 -m pele_platform.main {}\n'.format(
                          self.yaml)]
            slurm.writelines(lines2)


def create_20sbatch(ligchain, ligname, atom1, atom2, cpus=24, folder="pdb_files", test=False, initial=None,
                    file_list=None, cu=False):
    """
    creates for each of the mutants the yaml and slurm files
    ligchain (str): the chain ID where the ligand is located
    ligname (str): the residue name of the ligand in the PDB
    atom1 (str): atom of the residue to follow  --> chain ID:position:atom name
    atom2 (str): atom of the ligand to follow  --> chain ID:position:atom name
    cpus (str or int): how many cpus do you want to use
    """
    if not os.path.exists(folder):
        raise IOError("No directory named {}".format(folder))

    slurm_files = []
    if not file_list:
        for files in glob("{}/*.pdb".format(folder)):
            name = basename(files)
            name = name.replace(".pdb", "")
            run = CreateLaunchFiles(files, ligchain, ligname, atom1, atom2, cpus, test=test, initial=initial, cu=cu)
            run.match_dist()
            run.input_creation(name)
            run.slurm_creation(name)
            slurm_files.append(run.slurm)

    else:
        for files in file_list:
            name = basename(files)
            name = name.replace(".pdb", "")
            run = CreateLaunchFiles(files, ligchain, ligname, atom1, atom2, cpus, test=test, initial=initial, cu=cu)
            run.match_dist()
            run.input_creation(name)
            run.slurm_creation(name)
            slurm_files.append(run.slurm)

    return slurm_files


def main():
    folder, ligchain, ligname, atom1, atom2, cpus, test, cu = parse_args()
    slurm_files = create_20sbatch(ligchain, ligname, atom1, atom2, cpus=cpus, folder=folder, test=test, cu=cu)

    return slurm_files


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    yaml_list, slurm_list = main()
