"""
This script is designed to generate the sbatch file for the simulations to run.
"""

__author__ = "Ruite Xiang"
__license__ = "MIT"
__maintainer__ = "Ruite Xiang"
__email__ = "ruite.xiang@bsc.es"


import argparse
from subprocess import call
from os.path import basename
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
    parser.add_argument("--sbatch", required=False, action="store_false",
                        help="True if you want to lanch the simulation right after creating the slurm file")

    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atom1, args.atom2, args.cpus, args.test,
            args.cu, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen, args.consec, args.sbatch]


class CreateSlurmFiles:
    """
    Creates the 2 necessary files for the pele simulations
    """
    def __init__(self, input_, ligchain, ligname, atom1, atom2, length, position, cpus=24, dir_=None, hydrogen=True,
                  multiple=False, pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345, nord=False):
        """
        Initialize the CreateLaunchFiles object

        Parameters
        ___________
        input_: str
            A  PDB file's path
        ligchain: str
            the chain ID where the ligand is located
        ligname: str
            the residue name of the ligand in the PDB
        atom1: str
            atom of the residue to follow in this format --> chain ID:position:atom name
        atom2: str
            atom of the ligand to follow in this format --> chain ID:position:atom name
        length: int
            To calculate the real cpus necessary to run all the simulations
        position: list[str]
            [chain ID:position] of the residue, for example [A:139,..]
        cpus: int, optional
            how many cpus do you want to use
        dir_: str, optional
            Name of the folder ofr the simulations
        hydrogens: bool, optional
            Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
        multiple: bool, optional
            Specify if to mutate 2 positions at the same pdb
        pdb_dir: str, optional
            The name of the folder where the mutated PDB files will be stored
        consec: bool, optional
            Consecutively mutate the PDB file for several rounds
        test: bool, optional
            Setting the simulation to test mode
        cu: bool, optional
            Set it to true if there are coppers in the system
        seed: int, optional
            A seed number to make the simulations reproducible
        nord: bool, optional
            True if the system is managed by LSF
        """

        self.input = input_
        self.ligchain = ligchain
        self.ligname = ligname
        self.atom1 = atom1
        self.atom2 = atom2
        self.cpus = cpus
        self.test = test
        self.slurm = None
        self.cu = cu
        self.seed = seed
        self.nord = nord
        self.len = length
        self.position = position
        self.hydrogen = hydrogen
        self.multiple = multiple
        self.consec = consec
        self.dir = dir_
        self.pdb_dir = pdb_dir

    def slurm_creation(self):
        """
        Creates the slurm running files for PELE in sbatch managed systems
        """
        name = basename(self.input).replace(".pdb", "")
        self.slurm = "{}.sh".format(name)
        with open(self.slurm, "w") as slurm:
            lines = ["#!/bin/bash\n", "#SBATCH -J PELE\n", "#SBATCH --output={}.out\n".format(name),
                     "#SBATCH --error={}.err\n".format(name)]
            if self.test:
                lines.append("#SBATCH --qos=debug\n")
                self.cpus = 5
                real_cpus = self.cpus * self.len
                lines.append("#SBATCH --ntasks={}\n\n".format(real_cpus))
            else:
                real_cpus = self.cpus * self.len
                lines.append("#SBATCH --ntasks={}\n\n".format(real_cpus))

            lines2 = ['module purge\n',
                      'export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"\n',
                      'export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"\n',
                      'export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH\n',
                      'module load intel mkl impi gcc # 2> /dev/null\n', 'module load boost/1.64.0\n']

            argument_list = []
            posi = " ".join(self.position)
            arguments = "--input {} --position {} --ligchain {} --ligname {} --atom1 {} --atom2 {} ".format(
                self.input, posi, self.ligchain, self.ligname, self.atom1, self.atom2)
            argument_list.append(arguments)

            if self.seed != 12345:
                argument_list.append("--seed {} ".format(self.seed))
            if self.cpus != 24:
                argument_list.append("--cpus {} ".format(self.cpus))
            if self.hydrogen:
                argument_list.append("--hydrogen ")
            if self.consec:
                argument_list.append("--consec ")
            if self.multiple:
                argument_list.append("--multiple ")
            if self.cu:
                argument_list.append("--cu ")
            if self.nord:
                argument_list.append("--nord ")
            if self.pdb_dir != "pdb_files":
                argument_list.append("--pdb_dir {} ".format(self.pdb_dir))
            if self.dir:
                argument_list.append("--dir {} ".format(self.dir))
            if self.test:
                argument_list.append("--test ")

            all_arguments = "".join(argument_list)
            python = "/home/bsc72/bsc72661/.conda/envs/saturated/bin/python -m saturated_mutagenesis.simulation {}".format(
                all_arguments)
            lines2.append(python)
            lines.extend(lines2)
            slurm.writelines(lines)

        return self.slurm

    def slurm_nord(self, slurm_name):
        """
        Create slurm files for PELE in LSF managed systems

        Parameters
        ___________
        slurm_name: str
            Name of the file created
        """
        self.slurm = "{}.sh".format(slurm_name)
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
                      'python -m pele_platform.main {}\n']

            lines.extend(lines2)
            slurm.writelines(lines)


def main():
    input_, position, ligchain, ligname, atom1, atom2, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, \
    hydrogen, consec, sbatch = parse_args()

    if multiple and len(position) == 2:
        length = 400
    else:
        length = len(position) * 19 + 1
    run = CreateSlurmFiles(input_, ligchain, ligname, atom1, atom2, length, position, cpus, dir_, hydrogen,
                  multiple, pdb_dir, consec, test, cu, seed, nord)
    slurm = run.slurm_creation()
    if sbatch:
        call(["sbatch", "{}".format(slurm)])

if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
