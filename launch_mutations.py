import argparse
import os
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Make predictions")
    # main required arguments
    parser.add_argument("--folder", required=True, help="Include the folder where the pdb files are located")
    parser.add_argument("--chain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--resname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format: chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format: chain ID:position:atom name")
    #arguments = vars(parser.parse_args())
    args = parser.parse_args()
    return args.folder, args.chain, args.resname, args.atom1, args.atom2

class CreateLaunchFiles():
    def __init__(self, input_, chain, resname, atom1, atom2, cpus=24):
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
        self.yaml = None
        self.slurm = None

    def input_creation(self, yaml_name, test=False):
        """ create the .yaml input files for PELE"""
        if not os.path.exists("yaml_files"):
            os.mkdir("yaml_files")
        self.yaml = "input_files/{}.yaml".format(yaml_name)
        with open(self.yaml, "w") as inp:
            inp.write("system: '{}'".format(self.input))
            inp.write("chain: '{}'".format(self.chain))
            inp.write("resname: '{}'".format(self.resname))
            inp.write("induced_fit_exhaustive: true")
            inp.write("seed: 12345")
            inp.write("cpus: {}".format(self.cpus))
            inp.write("atom_dist:\n- '{}'\n- '{}'".format(self.atom1, self.atom2))
            inp.write("skip_preprocess: true")
            if test:
                inp.write("test: true")



    def slurm_creation(self, slurm_name, yaml_file, test=False):
        """Creates the slurm runing files for PELE"""
        if not os.path.exists("slurm_files"):
            os.mkdir("slurm_files")
        self.slurm = "slurm_files/{}.sh".format(slurm_name)
        with open(self.slurm, "w") as slurm:
            slurm.write("#!/bin/bash")
            slurm.write("#SBATCH -J PELE")
            slurm.write("#SBATCH --output=mpi_%j.out")
            slurm.write("#SBATCH --error=mpi_%j.err")
            slurm.write("#SBATCH --ntasks=50")
            if test:
                slurm.write("#SBATCH --qos=debug\n")
            slurm.write('module purgeexport PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"')
            slurm.write('export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"')
            slurm.write('export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH')
            slurm.write('module load impi')
            slurm.write('module load intel mkl impi gcc # 2> /dev/null')
            slurm.write('module load boost/1.64.0')
            slurm.write('/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8 -m pele_platform.main {}'.format(yaml_file))



def create_20sbatch(chain, resname, atom1, atom2, cpus=24, folder="pdb_files"):
    """
    creates for each of the mutants the yaml and slurm files

        chain: (str) the chain ID where the ligand is located
        resname: (str) the residue name of the ligand in the PDB
        atom1: (str) atom of the residue to follow in this format --> chain ID:position:atom name
        atom2: (str) atom of the ligand to follow in this format --> chain ID:position:atom name
        cpus: (str or int) how many cpus do you want to use

    """
    if not os.path.exists(folder):
        raise IOError("No directory named {}, run mutate_pdb.py first".format(folder))

    yaml_files = []
    slurm_files = []
    for file in glob.glob("{}/*.pdb".format(folder)):
        name = file.replace("{}/".format(folder), "")
        run = CreateLaunchFiles(file, chain, resname, atom1, atom2, cpus)
        run.input_creation(name)
        run.slurm_creation(name, yaml_file=run.yaml)
        yaml_files.append(run.yaml)
        slurm_files.append(run.slurm)

    return yaml_files, slurm_files

def main():
    folder, chain, resname, atom1, atom2 = parse_args()
    yaml_files, slurm_files = create_20sbatch(chain, resname, atom1, atom2, cpus=24, folder=folder)

    return yaml_files, slurm_files

if __name__ == "__main__":
    #Run this if this file is executed from command line but not if is imported as API
    yaml_files, slurm_files = main()
