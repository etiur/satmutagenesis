import argparse
import os
import glob
from helper import map_atom_string

def parse_args():
    parser = argparse.ArgumentParser(description="Generate running files for PELE")
    # main required arguments
    parser.add_argument("--folder", required=False, default="pdb_files",
                        help="Include the folder where the pdb files are located")
    parser.add_argument("--chain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--resname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=24, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--cu", required=False, action="store_true")
    parser.add_argument("--test", required=False, action="store_true")
    args = parser.parse_args()
    return args.folder, args.chain, args.resname, args.atom1, args.atom2, args.cpus, args.test, args.cu

class CreateLaunchFiles():
    def __init__(self, input_, chain, resname, atom1, atom2, cpus=24, test=False, initial=None, cu=False):
        """
        input_: (str) PDB files path
        chain: (str) the chain ID where the ligand is located
        resname: (str) the residue name of the ligand in the PDB
        atom1: (str) atom of the residue to follow in this format --> chain ID:position:atom name
        atom2: (str) atom of the ligand to follow in this format --> chain ID:position:atom name
        cpus: (str or int) how many cpus do you want to use
        """
        self.input = input_
        self.chain = chain
        self.resname = resname
        self.atom1 = atom1
        self.atom2 = atom2
        self.cpus = cpus
        self.test = test
        self.yaml = None
        self.slurm = None
        self.initial = initial
        self.cu = cu

    def match_dist(self):
        """ match the user coordinates to pmx PDB coordinates"""
        if self.initial:
            self.atom1 = map_atom_string(self.atom1, self.initial, self.input)
            self.atom2 = map_atom_string(self.atom2, self.initial, self.input)
        else:
            pass

    def input_creation(self, yaml_name):
        """ create the .yaml input files for PELE"""
        if not os.path.exists("yaml_files"):
            os.mkdir("yaml_files")
        self.yaml = "yaml_files/{}.yaml".format(yaml_name)
        with open(self.yaml, "w") as inp:
            inp.write("system: '{}'\n".format(self.input))
            inp.write("chain: '{}'\n".format(self.chain))
            inp.write("resname: '{}'\n".format(self.resname))
            inp.write("induced_fit_exhaustive: true\n")
            inp.write("seed: 12345\n")
            if self.test:
                inp.write("test: true\n")
                self.cpus = 5
            inp.write("cpus: {}\n".format(self.cpus))
            inp.write("atom_dist:\n- '{}'\n- '{}'\n".format(self.atom1, self.atom2))
            if self.cu:
                path = "/gpfs/projects/bsc72/ruite/examples/cuz"
                inp.write("templates:\n- '{}'\n".format(path))
            inp.write("pele_license: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/license'\n")
            inp.write("pele_exec: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/bin/PELE-1.6.1_mpi'\n")

    def slurm_creation(self, slurm_name):
        """Creates the slurm runing files for PELE"""
        if not os.path.exists("slurm_files"):
            os.mkdir("slurm_files")
        self.slurm = "slurm_files/{}.sh".format(slurm_name)
        with open(self.slurm, "w") as slurm:
            slurm.write("#!/bin/bash\n")
            slurm.write("#SBATCH -J PELE\n")
            slurm.write("#SBATCH --output={}.out\n".format(slurm_name))
            slurm.write("#SBATCH --error={}.err\n".format(slurm_name))
            if self.test:
                slurm.write("#SBATCH --qos=debug\n")
                self.cpus = 5
            slurm.write("#SBATCH --ntasks={}\n\n".format(self.cpus))
            slurm.write('module purge\n')
            slurm.write('export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"\n')
            slurm.write('export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"\n')
            slurm.write('export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH\n')
            slurm.write('module load impi\n')
            slurm.write('module load intel mkl impi gcc # 2> /dev/null\n')
            slurm.write('module load boost/1.64.0\n')
            slurm.write('/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8 -m pele_platform.main {}\n'.format(self.yaml))

def create_20sbatch(chain, resname, atom1, atom2, cpus=24, folder="pdb_files", test=False, initial=None, file_list=None, cu=False):
    """
    creates for each of the mutants the yaml and slurm files

        chain: (str) the chain ID where the ligand is located
        resname: (str) the residue name of the ligand in the PDB
        atom1: (str) atom of the residue to follow  --> chain ID:position:atom name
        atom2: (str) atom of the ligand to follow  --> chain ID:position:atom name
        cpus: (str or int) how many cpus do you want to use
    """
    if not os.path.exists(folder):
        raise IOError("No directory named {}".format(folder))

    slurm_files = []
    if not file_list:
        for file in glob.glob("{}/*.pdb".format(folder)):
            name = file.replace("{}/".format(folder), "")
            name = name.replace("{}".format(folder), "")
            name = name.replace(".pdb", "")
            run = CreateLaunchFiles(file, chain, resname, atom1, atom2, cpus, test=test, initial=initial, cu=cu)
            run.match_dist()
            run.input_creation(name)
            run.slurm_creation(name)
            slurm_files.append(run.slurm)
    else:
        for file in file_list:
            name = file.replace("{}/".format(file.split("/")[0]), "")
            name = name.replace(".pdb", "")
            run = CreateLaunchFiles(file, chain, resname, atom1, atom2, cpus, test=test, initial=initial, cu=cu)
            run.match_dist()
            run.input_creation(name)
            run.slurm_creation(name)
            slurm_files.append(run.slurm)

    return slurm_files

def main():
    folder, chain, resname, atom1, atom2, cpus, test, cu = parse_args()
    yaml_files, slurm_files = create_20sbatch(chain, resname, atom1, atom2, cpus=cpus, folder=folder, test=test, cu=cu)

    return yaml_files, slurm_files

if __name__ == "__main__":
    #Run this if this file is executed from command line but not if is imported as API
    yaml_files, slurm_files = main()
