import argparse
import os
from helper import map_atom_string, isiterable
from os.path import basename, join, isfile, isdir


def parse_args():
    parser = argparse.ArgumentParser(description="Generate running files for PELE")
    # main required arguments
    parser.add_argument("--folder", required=True,
                        help="Include the folder where the pdb files are located")
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
    parser.add_argument("--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    args = parser.parse_args()

    return [args.folder, args.ligchain, args.ligname, args.atom1, args.atom2, args.cpus, args.test, args.cu,
            args.seed, args.nord]


class CreateLaunchFiles:
    def __init__(self, input_, ligchain, ligname, atom1, atom2, cpus=24,
                 test=False, initial=None, cu=False, seed=12345, nord=False):
        """
        input_ (str): PDB files path
        ligchain (str): the chain ID where the ligand is located
        ligname (str): the residue name of the ligand in the PDB
        atom1 (str): atom of the residue to follow in this format --> chain ID:position:atom name
        atom2 (str): atom of the ligand to follow in this format --> chain ID:position:atom name
        cpus (str or int): How many cpus do you want to use
        test (boolean): Setting the simulation to test mode
        initial (file): The initial PDB file before the modification by pmx
        cu (boolean): Set it to true if there are coppers in the system
        seed (int): A seed number to make the simulations reproducible
        nord (boolean): True if the system is managed by LSF
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
        self.seed = seed
        self.nord = nord

    def _match_dist(self):
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
        yaml_name (str): Name for the input file for the simulation
        """
        self._match_dist()
        if not os.path.exists("yaml_files"):
            os.mkdir("yaml_files")
        self.yaml = "yaml_files/{}.yaml".format(yaml_name)
        with open(self.yaml, "w") as inp:
            lines = ["system: '{}'\n".format(self.input), "chain: '{}'\n".format(self.ligchain),
                     "resname: '{}'\n".format(self.ligname), "induced_fit_exhaustive: true\n",
                     "seed: {}\n".format(self.seed)]
            if not self.nord:
                lines.append("usesrun: true\n")
            if yaml_name != "original":
                lines.append("working_folder: {}/PELE_{}\n".format(yaml_name[:-1], yaml_name))
            else:
                lines.append("working_folder: PELE_{}\n".format(yaml_name))
            if self.test:
                lines.append("test: true\n")
                self.cpus = 5
            lines2 = ["cpus: {}\n".format(self.cpus), "atom_dist:\n- '{}'\n- '{}'\n".format(self.atom1, self.atom2),
                      "pele_license: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/license'\n",
                      "pele_exec: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/bin/PELE-1.6.1_mpi'\n"]
            if self.cu:
                path = "/gpfs/projects/bsc72/ruite/examples/cuz"
                lines2.append("templates:\n- '{}'\n".format(path))
            lines.extend(lines2)
            inp.writelines(lines)

    def slurm_creation(self, slurm_name):
        """
        Creates the slurm running files for PELE in sbatch managed systems
        slurm_name (str): Name for the batch file
        """
        if not os.path.exists("slurm_files"):
            os.mkdir("slurm_files")
        self.slurm = "slurm_files/{}.sh".format(slurm_name)
        with open(self.slurm, "w") as slurm:
            lines = ["#!/bin/bash\n", "#SBATCH -J PELE\n", "#SBATCH --output={}.out\n".format(slurm_name),
                     "#SBATCH --error={}.err\n".format(slurm_name)]
            if self.test:
                lines.append("#SBATCH --qos=debug\n")
                self.cpus = 5
                lines.append("#SBATCH --ntasks={}\n\n".format(self.cpus))
            else:
                lines.append("#SBATCH --ntasks={}\n\n".format(self.cpus))

            lines2 = ['module purge\n',
                      'export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"\n',
                      'export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"\n',
                      'export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH\n',
                      'module load intel mkl impi gcc # 2> /dev/null\n', 'module load boost/1.64.0\n',
                      '/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8 -m pele_platform.main {}\n'.format(
                          self.yaml)]

            lines.extend(lines2)
            slurm.writelines(lines)

    def slurm_nord(self, slurm_name):
        """
        Create slurm files for PELE in LSF managed systems
        slurm_name (str): Name of the file created
        """
        if not os.path.exists("slurm_files"):
            os.mkdir("slurm_files")
        self.slurm = "slurm_files/{}.sh".format(slurm_name)
        with open(self.slurm, "w") as slurm:
            lines = ["#!/bin/bash\n", "#BSUB -J PELE\n", "#BSUB -oo {}.out\n".format(slurm_name),
                     "#BSUB -eo {}.err\n".format(slurm_name)]
            if self.test:
                lines.append("#BSUB -q debug\n")
                self.cpus = 5
                lines.append("#BSUB -W 01:00\n")
                lines.append("#BSUB -n {}\n\n".format(self.cpus))
            else:
                lines.append("#BSUB -W 48:00\n")
                lines.append("#BSUB -n {}\n\n".format(self.cpus))

            lines2 = ['module purge\n',
                      'module load intel gcc/latest openmpi/1.8.1 boost/1.63.0 PYTHON/3.7.4 MKL/11.3 GTK+3/3.2.4\n',
                      'export PYTHONPATH=/gpfs/projects/bsc72/PELEPlatform/1.5.1/pele_platform:$PYTHONPATH\n',
                      'export PYTHONPATH=/gpfs/projects/bsc72/PELEPlatform/1.5.1/dependencies:$PYTHONPATH\n',
                      'export PYTHONPATH=/gpfs/projects/bsc72/adaptiveSampling/bin_nord/v1.6.2/:$PYTHONPATH\n',
                      'export PYTHONPATH=/gpfs/projects/bsc72/PELEPlatform/external_deps/:$PYTHONPATH\n',
                      'export PYTHONPATH=/gpfs/projects/bsc72/lib/site-packages_mn3:$PYTHONPATH\n',
                      'export MPLBACKEND=Agg\n', 'export OMPI_MCA_coll_hcoll_enable=0\n', 'export OMPI_MCA_mtl=^mxm\n'
                      'python -m pele_platform.main {}\n'.format(self.yaml)]

            lines.extend(lines2)
            slurm.writelines(lines)


def create_20sbatch(ligchain, ligname, atom1, atom2, file_, cpus=24, test=False, initial=None,
                    cu=False, seed=12345, nord=False):
    """
    creates for each of the mutants the yaml and slurm files
    ligchain (str): the chain ID where the ligand is located
    ligname (str): the residue name of the ligand in the PDB
    atom1 (str): atom of the residue to follow  --> chain ID:position:atom name
    atom2 (str): atom of the ligand to follow  --> chain ID:position:atom name
    cpus (str or int): how many cpus do you want to use
    test (boolean): Setting the simulation to test mode
    initial (file): The initial PDB file before the modification by pmx
    file_ (iterable, folder or a file): An iterable of the location of the different pdb files, a name of the folder
    or a file of the path to the different pdb files
    cu (boolean): Set it to true if there are coppers in the system
    seed (int): A seed number to make the simulations reproducible
    nord (boolean): True if the system is managed by LSF
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
    for files in file_list:
        files = files.strip("\n")
        name = basename(files)
        name = name.replace(".pdb", "")
        run = CreateLaunchFiles(files, ligchain, ligname, atom1, atom2, cpus, test=test,
                                initial=initial, cu=cu, seed=seed, nord=nord)
        run.input_creation(name)
        if not nord:
            run.slurm_creation(name)
        else:
            run.slurm_nord(name)
        slurm_files.append(run.slurm)

    return slurm_files


def main():
    folder, ligchain, ligname, atom1, atom2, cpus, test, cu, seed, nord = parse_args()
    slurm_files = create_20sbatch(ligchain, ligname, atom1, atom2,
                                  cpus=cpus, file_=folder, test=test, cu=cu, seed=seed, nord=nord)

    return slurm_files


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    yaml_list, slurm_list = main()
