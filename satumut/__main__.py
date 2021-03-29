"""
This script is designed to generate the sbatch file for the simulations to run.
"""

__author__ = "Ruite Xiang"
__license__ = "MIT"
__maintainer__ = "Ruite Xiang"
__email__ = "ruite.xiang@bsc.es"

import argparse
import os
from os.path import basename
from helper import neighbourresidues
from Bio import PDB


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("-i", "--input", required=True, nargs="+", help="Include one or more PDB path, all the"
                                                                        "rest remains the same for both PDBs")
    parser.add_argument("-p", "--position", required=False, nargs="+",
                        help="Include one or more chain IDs and positions -> Chain ID:position")
    parser.add_argument("-lc", "--ligchain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("-ln", "--ligname", required=True, help="The ligand residue name")
    parser.add_argument("-at", "--atoms", required=True, nargs="+",
                        help="Series of atoms of the residues to follow in this format -> chain ID:position:atom name")
    parser.add_argument("-cpm", "--cpus_per_mutant", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-tcpus", "--total_cpus", required=False, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-po", "--polarize_metals", required=False, action="store_true",
                        help="used if there are metals in the system")
    parser.add_argument("-fa", "--polarization_factor", required=False, type=int,
                        help="The number to divide the charges")
    parser.add_argument("-t", "--test", required=False, action="store_true",
                        help="Used if you want to run a test before")
    parser.add_argument("-n", "--nord", required=False, action="store_true",
                        help="used if LSF is the utility managing the jobs")
    parser.add_argument("-m", "--multiple", required=False, action="store_true",
                        help="if you want to mutate 2 residue in the same pdb")
    parser.add_argument("-s", "--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    parser.add_argument("-d", "--dir", required=False,
                        help="The name of the folder for all the simulations")
    parser.add_argument("-pd", "--pdb_dir", required=False, default="pdb_files",
                        help="The name for the mutated pdb folder")
    parser.add_argument("-hy", "--hydrogen", required=False, action="store_false", help="leave it to default")
    parser.add_argument("-co", "--consec", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")
    parser.add_argument("-sb", "--sbatch", required=False, action="store_false",
                        help="True if you want to lanch the simulation right after creating the slurm file")
    parser.add_argument("-st", "--steps", required=False, type=int, default=800,
                        help="The number of PELE steps")
    parser.add_argument("--dpi", required=False, default=800, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("--box", required=False, default=30, type=int,
                        help="Set how many data points are used for the boxplot")
    parser.add_argument("-tr", "--trajectory", required=False, default=10, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--out", required=False, default="summary",
                        help="Name of the summary file created at the end of the analysis")
    parser.add_argument("--plot", required=False,
                        help="Path of the plots folder")
    parser.add_argument("-an", "--analyse", required=False, choices=("energy", "distance", "both"), default="distance",
                        help="The metric to measure the improvement of the system")
    parser.add_argument("--thres", required=False, default=-0.1, type=float,
                        help="The threshold for the improvement which will affect what will be included in the summary")
    parser.add_argument("-sm", "--single_mutagenesis", required=False,
                        help="Specifiy the name of the residue that you want the "
                             "original residue to be mutated to. Both 3 letter "
                             "code and 1 letter code can be used. You can even specify the protonated states")
    parser.add_argument("-PR", "--plurizyme_at_and_res", required=False,
                        help="Specify the chain ID, residue number and the PDB atom name that"
                             "will set the list of the neighbouring residues for the"
                             "next round. Example: chain ID:position:atom name")
    parser.add_argument("-r", "--radius", required=False, default=5.0, type=float,
                        help="The radius around the selected atom to search for the other residues")
    parser.add_argument("-f", "--fixed_resids", required=False, default=(), nargs='+',
                        help="Specify the list of residues that you don't want"
                             "to have mutated (Must write the list of residue position"
                             "numbers)")
    parser.add_argument("-cpt", "--cpus_per_task", required=False, default=1, type=int,
                        help="Include the number of cpus per task desired")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atoms, args.cpus_per_mutant, args.test,
            args.polarize_metals, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen,
            args.consec, args.sbatch, args.steps, args.dpi, args.box, args.trajectory, args.out, args.plot, args.analyse,
            args.thres, args.single_mutagenesis, args.plurizyme_at_and_res, args.radius, args.fixed_resids,
            args.cpus_per_task, args.polarization_factor, args.total_cpus, args.xtc, args.catalytic_distance]


class CreateSlurmFiles:
    """
    Creates the 2 necessary files for the pele simulations
    """

    def __init__(self, input_, ligchain, ligname, atoms, position=(), cpus_mutant=25, dir_=None, hydrogen=True,
                 multiple=False, pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345, nord=False,
                 steps=800, dpi=800, box=30, traj=10, output="summary", plot_dir=None, opt="distance", thres=-0.1,
                 single_mutagenesis=None, plurizyme_at_and_res=None, radius=5.0, fixed_resids=(), cpus_task=1,
                 factor=None, total_cpus=None, xtc=False, cata_dist=3.5):
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
        atoms: list[str]
            list of atom of the residue to follow, in this format --> chain ID:position:atom name
        position: list[str]
            [chain ID:position] of the residue, for example [A:139,..]
        cpus_mutant: int, optional
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
            Set it to true if there are metals in the system
        seed: int, optional
            A seed number to make the simulations reproducible
        nord: bool, optional
            True if the system is managed by LSF
        steps: int, optional
            The number of PELE steps
        dpi : int, optional
            The quality of the plots
        box : int, optional
            how many points are used for the box plots
        traj : int, optional
            how many top pdbs are extracted from the trajectories
        output : str, optional
            name of the output file for the pdfs
        plot_dir : str
            Name for the results folder
        opt : str, optional
            choose if to analyse distance, energy or both
        thres : float, optional
            The threshold for the mutations to be included in the pdf
        single_mutagenesis: str
            The new residue to mutate the positions to, in 3 letter or 1 letter code
        plurizyme_at_and_res: str
            Chain_ID:position:atom_name, which will be the center around the search for neighbours
        radius: float, optional
            The radius for the neighbours search
        fixed_resids. list[position_num]
            A list of residues positions to avoid mutating
        factor: int, optional
            The number to divide the metal charges
        Total_cpus: int, optional
            The total number of cpus available
        xtc: bool, optional
            Set to True if you want to change the pdb format to xtc
        cata_dist: float, optional
            The catalytic distance
        """
        assert len(atoms) % 2 == 0, "Introduce pairs of atoms to follow"
        self.input = input_
        self.ligchain = ligchain
        self.ligname = ligname
        self.atoms = " ".join(atoms)
        self.cpus = cpus_mutant
        self.test = test
        self.slurm = None
        self.cu = cu
        self.seed = seed
        self.nord = nord
        if multiple and len(position) == 2:
            self.len = 400
        elif single_mutagenesis and plurizyme_at_and_res:
            _ = neighbourresidues(input_, plurizyme_at_and_res, radius, fixed_resids)
            self.len = len(_)
        else:
            self.len = len(position) * 19 + 1
        if len(position) != 0:
            self.position = " ".join(position)
        else:
            self.position = None
        self.hydrogen = hydrogen
        self.multiple = multiple
        self.consec = consec
        self.dir = dir_
        self.pdb_dir = pdb_dir
        self.steps = steps
        self.dpi = dpi
        self.box = box
        self.traj = traj
        self.output = output
        self.plot_dir = plot_dir
        self.opt = opt
        self.thres = thres
        self.single = single_mutagenesis
        self.pluri = plurizyme_at_and_res
        self.radius = radius
        self.avoid = fixed_resids
        self.cpus_task = cpus_task
        self.factor = factor
        self.total_cpus = total_cpus
        self.xtc = xtc
        self.cata_dist = cata_dist

    def _size(self):
        """
        A function to calculate the size of the PDB file

        Returns
        ________
        residues length: int
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(self.input[:-4], self.input)
        residues = list(structure.get_residues())

        return len(residues)

    def slurm_creation(self):
        """
        Creates the slurm running files for PELE in sbatch managed systems
        """
        if not self.dir:
            name = basename(self.input).replace(".pdb", "")
        else:
            name = basename(self.dir)
        self.slurm = "{}.sh".format(name)
        with open(self.slurm, "w") as slurm:
            lines = ["#!/bin/bash\n", "#SBATCH -J {}\n".format(name), "#SBATCH --output={}.out\n".format(name),
                     "#SBATCH --error={}.err\n".format(name)]
            if self.test:
                lines.append("#SBATCH --qos=debug\n")
                self.cpus = 5
                real_cpus = self.cpus * self.len
                lines.append("#SBATCH --ntasks={}\n\n".format(real_cpus))
            else:
                if self.total_cpus:
                    real_cpus = self.total_cpus
                else:
                    real_cpus = self.cpus * self.len + 2
                lines.append("#SBATCH --ntasks={}\n".format(real_cpus))
                lines.append("#SBATCH --cpus-per-task={}\n\n".format(self.cpus_task))

            lines2 = ['module purge\n',
                      'export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"\n',
                      'export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"\n',
                      'module load intel mkl impi gcc # 2> /dev/null\n', 'module load boost/1.64.0 ANACONDA/2019.10\n',
                      'eval "$(conda shell.bash hook)"\n',
                      "conda activate /gpfs/projects/bsc72/conda_envs/platform/1.6.0b\n\n"]

            argument_list = []
            arguments = "-i {} -lc {} -ln {} -at {} ".format(self.input, self.ligchain, self.ligname, self.atoms)
            argument_list.append(arguments)
            if self.position:
                argument_list.append("-p {} ".format(self.position))
            if self.seed != 12345:
                argument_list.append("--seed {} ".format(self.seed))
            if self.cpus != 25:
                argument_list.append("-cpm {} ".format(self.cpus))
            if self.total_cpus:
                argument_list.append("-tcpus {} ".format(self.total_cpus))
            if not self.hydrogen:
                argument_list.append("-hy ")
            if self.consec:
                argument_list.append("-co ")
            if self.multiple:
                argument_list.append("-m ")
            if self.cu:
                argument_list.append("-po ")
            if self.nord:
                argument_list.append("--nord ")
            if self.pdb_dir != "pdb_files":
                argument_list.append("-pd {} ".format(self.pdb_dir))
            if self.dir:
                argument_list.append("--dir {} ".format(self.dir))
            if self.test:
                argument_list.append("--test ")
            if self.xtc:
                argument_list.append("-x ")
            if self.steps != 800:
                argument_list.append("--steps {} ".format(self.steps))
            if self.dpi != 800:
                argument_list.append("--dpi {} ".format(self.dpi))
            if self.box != 30:
                argument_list.append("--box {} ".format(self.box))
            if self.traj != 10:
                argument_list.append("-tr {} ".format(self.traj))
            if self.output != "summary":
                argument_list.append("--out {} ".format(self.output))
            if self.plot_dir:
                argument_list.append("--plot {} ".format(self.plot_dir))
            if self.opt != "distance":
                argument_list.append("-an {} ".format(self.opt))
            if self.thres != -0.1:
                argument_list.append("--thres {} ".format(self.thres))
            if self.cata_dist != 3.5:
                argument_list.append("-cd {} ".format(self.cata_dist))
            if self.single and self.pluri:
                argument_list.append("-sm {} ".format(self.single))
                argument_list.append("-PR {} ".format(self.pluri))
                if self.radius != 5.0:
                    argument_list.append("-r {} ".format(self.radius))
                if len(self.avoid) != 0:
                    argument_list.append("-f {} ".format(self.avoid))
            if self.cu and self.factor:
                argument_list.append("-fa {} ".format(self.factor))
            all_arguments = "".join(argument_list)
            python = "/gpfs/projects/bsc72/conda_envs/saturated/bin/python -m satumut.simulation {}\n".format(
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
    input_, position, ligchain, ligname, atoms, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, \
    hydrogen, consec, sbatch, steps, dpi, box, traj, out, plot_dir, analysis, thres, single_mutagenesis, \
    plurizyme_at_and_res, radius, fixed_resids, cpus_per_task, factor, total_cpus, xtc, cata_dist = parse_args()
    if dir_ and len(input_) > 1:
        dir_ = None
    for inp in input_:
        run = CreateSlurmFiles(inp, ligchain, ligname, atoms, position, cpus, dir_, hydrogen,
                               multiple, pdb_dir, consec, test, cu, seed, nord, steps, dpi, box, traj,
                               out, plot_dir, analysis, thres, single_mutagenesis, plurizyme_at_and_res, radius,
                               fixed_resids, cpus_per_task, factor, total_cpus, xtc)
        slurm = run.slurm_creation()
        if sbatch:
            os.system("sbatch {}".format(slurm))


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
