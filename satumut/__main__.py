"""
This script is designed to generate the sbatch file for the simulations to run.
"""

__author__ = "Ruite Xiang"
__license__ = "MIT"
__maintainer__ = "Ruite Xiang"
__email__ = "ruite.xiang@bsc.es"

import argparse
from subprocess import call
from pathlib import Path
from .helper import neighbourresidues
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
    parser.add_argument("-at", "--atoms", required=False, nargs="+",
                        help="Series of atoms of the residues to follow in this format -> chain ID:position:atom name")
    parser.add_argument("-cpm", "--cpus_per_mutant", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-tcpus", "--total_cpus", required=False, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-cpt", "--cpus_per_task", required=False, type=int,
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
                        help="True if you want to launch the simulation right after creating the slurm file")
    parser.add_argument("-st", "--steps", required=False, type=int, default=1000,
                        help="The number of PELE steps")
    parser.add_argument("--dpi", required=False, default=800, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("-tr", "--trajectory", required=False, default=5, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--plot", required=False,
                        help="Path of the plots folder")
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
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    parser.add_argument("-tem", "--template", required=False, nargs="+",
                        help="Path to external forcefield templates")
    parser.add_argument("-sk", "--skip", required=False, nargs="+",
                        help="skip the processing of ligands by PlopRotTemp")
    parser.add_argument("-rot", "--rotamers", required=False, nargs="+",
                        help="Path to external rotamers templates")
    parser.add_argument("-e", "--equilibration", required=False, action="store_true",
                        help="Set equilibration")
    parser.add_argument("-l", "--log", required=False, action="store_true",
                        help="write logs")
    parser.add_argument("-da", "--dihedral_atoms", required=False, nargs="+",
                        help="The 4 atom necessary to calculate the dihedrals in format chain id:res number:atom name")
    parser.add_argument("-im", "--improve", required=False, choices=("R", "S"), default="R",
                        help="The enantiomer that should improve")
    parser.add_argument("-tu", "--turn", required=False, type=int,
                        help="the round of plurizyme generation, not needed for the 1st round")
    parser.add_argument("-en", "--energy_threshold", required=False, type=int,
                        help="An energy threshold that limits the points of scatter plots")
    parser.add_argument("--QM", required=False, help="The path to the QM charges")
    parser.add_argument("-br","--box_radius", required=False, type=int, help="Radius of the exploration box")
    parser.add_argument("-mut", "--mutation", required=False, nargs="+",
                        choices=('ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN',
                                 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR'),
                        help="The aminoacids to mutate to in 3 letter code")
    parser.add_argument("-cst", "--conservative", required=False, choices=(1, 2), default=None, type=int,
                        help="How conservative should the mutations be, choices are 1 (most conservative) and 2")
    parser.add_argument("-pw", "--profile_with", required=False, choices=("Binding Energy", "currentEnergy"),
                        default="Binding Energy", help="The metric to generate the pele profiles with")
    parser.add_argument("-w", "--wild", required=False, default=None,
                        help="The path to the folder where the reports from wild type simulation are")
    parser.add_argument("-scr", "--side_chain_resolution", required=False, type=int, default=10,
                        help="Affects the side chain sampling, the smaller the more accurate")
    parser.add_argument("-ep", "--epochs", required=False, type=int, default=1,
                        help="the number of adaptive epochs to run")
    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atoms, args.cpus_per_mutant, args.test,
            args.polarize_metals, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen,
            args.consec, args.sbatch, args.steps, args.dpi, args.trajectory, args.plot, args.single_mutagenesis,
            args.plurizyme_at_and_res, args.radius, args.fixed_resids, args.polarization_factor, args.total_cpus,
            args.xtc, args.catalytic_distance, args.template, args.skip, args.rotamers, args.equilibration, args.log,
            args.cpus_per_task, args.improve, args.turn, args.energy_threshold, args.QM, args.dihedral_atoms,
            args.box_radius, args.mutation, args.conservative, args.profile_with, args.wild, args.side_chain_resolution,
            args.epochs]


class CreateSlurmFiles:
    """
    Creates the 2 necessary files for the pele simulations
    """

    def __init__(self, input_, ligchain, ligname, atoms, position=(), cpus_mutant=25, dir_=None, hydrogen=True,
                 multiple=False, pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345, nord=False,
                 steps=1000, dpi=800, traj=5, plot_dir=None,  single_mutagenesis=None, plurizyme_at_and_res=None,
                 radius=5.0, fixed_resids=(), factor=None, total_cpus=None, xtc=False, cata_dist=3.5, template=None,
                 skip=None, rotamers=None, equilibration=True, log=False, cpt=None, improve="R", turn=None,
                 energy_thres=None, QM=None, dihedral=None, box_radius=None, mut=None, conservative=None,
                 profile_with="Binding Energy", wild=None, side_chain_resolution=10, epochs=1):
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
        traj : int, optional
            how many top pdbs are extracted from the trajectories
        plot_dir : str
            Name for the results folder
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
        template: str, optional
            Path to the external forcefield templates
        skip: str, optional
            Skip the processing of ligands by PlopRotTemp
        rotamers: str: optional
            Path to the external rotamers
        equilibration: bool, optional
            True to include equilibration steps
        log: bool, optional
            True to write the log files of simulations
        improve: str
            The enantiomer that should improve
        turn: int, optional
            The round of the plurizyme generation
        energy_thres: int, optional
            The binding energy to consider for catalytic poses
        QM: str, optional
            The path to the QM charges
        dihedral: list[str]
            The 4 atoms that form the dihedral in format chain ID:position:atom name
        box_radius: int, optional
            The radius of the exploration box
        mut: list[str], optional
            The list of mutations to perform
        conservative: int, optional
            How conservative should be the mutations according to Blossum62
        profile_with: str, optional
            The metric to generate the pele profiles with
        side_chain_resolution: int, optional
            The resolution of the side chain sampling, the smaller the better
        epochs: int, optional
            The number of adaptive epochs to run
        """

        self.input = Path(input_)
        self.ligchain = ligchain
        self.ligname = ligname
        if atoms:
            assert len(atoms) % 2 == 0, "The number of atoms to follow should be multiple of 2"
            self.atoms = " ".join(atoms)
        else:
            self.atoms = None
        self.cpus = cpus_mutant
        self.test = test
        self.slurm = None
        self.cu = cu
        self.seed = seed
        self.nord = nord
        if single_mutagenesis and plurizyme_at_and_res:
            _ = neighbourresidues(input_, plurizyme_at_and_res, radius, fixed_resids)
            self.len = len(_)
        else:
            self.len = 4
        if position:
            self.position = " ".join(position)
        else:
            self.position = None
        self.hydrogen = hydrogen
        self.multiple = multiple
        self.consec = consec
        self.dir = Path(dir_)
        self.pdb_dir = pdb_dir
        self.steps = steps
        self.dpi = dpi
        self.traj = traj
        self.plot_dir = plot_dir
        self.single = single_mutagenesis
        self.pluri = plurizyme_at_and_res
        self.radius = radius
        self.avoid = fixed_resids
        self.factor = factor
        self.total_cpus = total_cpus
        self.xtc = xtc
        self.cata_dist = cata_dist
        if template:
            self.template = " ".join(template)
        else:
            self.template = None
        if skip:
            self.skip = " ".join(skip)
        else:
            self.skip = None
        if rotamers:
            self.rotamer = " ".join(rotamers)
        else:
            self.rotamer = None
        self.equilibration = equilibration
        self.log = log
        if dihedral:
            self.dihedral_atoms = " ".join(dihedral)
        else:
            self.dihedral_atoms = None
        self.cpt = cpt
        self.improve = improve
        self.turn = turn
        self.energy_thres = energy_thres
        self.qm = QM
        self.box_radius = box_radius
        if mut:
            self.mut = " ".join(mut)
        else:
            self.mut = None
        self.conservative = conservative
        self.profile_with = profile_with
        self.wild = wild
        self.resolution = side_chain_resolution
        self.epochs = epochs

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
            name = self.input.stem
        else:
            name = self.dir.name
        self.slurm = f"{name}.sh"
        with open(self.slurm, "w") as slurm:
            lines = ["#!/bin/bash\n", f"#SBATCH -J {name}\n", f"#SBATCH --output={name}.out\n",
                     f"#SBATCH --error={name}.err\n"]
            if self.test:
                lines.append("#SBATCH --qos=debug\n")
            if self.cpt:
                lines.append(f"#SBATCH --cpus-per-task={self.cpt}\n")
            if self.total_cpus:
                real_cpus = self.total_cpus
            else:
                real_cpus = self.cpus * self.len + 1
            lines.append(f"#SBATCH --ntasks={real_cpus}\n\n")

            lines2 = ['module purge\n',
                      'module load intel mkl impi gcc # 2> /dev/null\n', 'module load boost/1.64.0 ANACONDA/2019.10\n',
                      'eval "$(conda shell.bash hook)"\n',
                      "conda activate /gpfs/projects/bsc72/conda_envs/platform/1.6.2\n\n"]

            argument_list = []
            arguments = f"-i {self.input} -lc {self.ligchain} -ln {self.ligname} "
            argument_list.append(arguments)
            if self.atoms:
                argument_list.append(f"-at {self.atoms} ")
            if self.position:
                argument_list.append(f"-p {self.position} ")
            if self.seed != 12345:
                argument_list.append(f"--seed {self.seed} ")
            if self.cpus != 25:
                argument_list.append(f"-cpm {self.cpus} ")
            if self.total_cpus:
                argument_list.append(f"-tcpus {self.total_cpus} ")
            if not self.hydrogen:
                argument_list.append("-hy ")
            if self.consec:
                argument_list.append("-co ")
            if self.multiple:
                argument_list.append("-m ")
            if self.cu:
                argument_list.append("-po ")
            if self.qm:
                argument_list.append(f"--QM {self.qm} ")
            if self.box_radius:
                argument_list.append(f"-br {self.box_radius} ")
            if self.nord:
                argument_list.append("--nord ")
            if self.pdb_dir != "pdb_files":
                argument_list.append(f"-pd {self.pdb_dir} ")
            if self.epochs != 1:
                argument_list.append(f"-ep {self.epochs} ")
            if self.dir:
                argument_list.append(f"--dir {self.dir} ")
            if self.equilibration:
                argument_list.append("-e ")
            if self.profile_with != "Binding Energy":
                argument_list.append(f"-pw {self.profile_with} ")
            if self.log:
                argument_list.append("-l ")
            if self.wild:
                argument_list.append(f"-w {self.wild} ")
            if self.xtc:
                argument_list.append("-x ")
            if self.steps != 1000:
                argument_list.append(f"--steps {self.steps} ")
            if self.dpi != 800:
                argument_list.append(f"--dpi {self.dpi} ")
            if self.traj != 5:
                argument_list.append(f"-tr {self.traj} ")
            if self.plot_dir:
                argument_list.append(f"--plot {self.plot_dir} ")
            if self.resolution:
                argument_list.append(f"-scr {self.resolution} ")
            if self.cata_dist != 3.5:
                argument_list.append(f"-cd {self.cata_dist} ")
            if self.single and self.pluri:
                argument_list.append(f"-sm {self.single} ")
                argument_list.append(f"-PR {self.pluri} ")
                if self.radius != 5.0:
                    argument_list.append(f"-r {self.radius} ")
                if len(self.avoid) != 0:
                    argument_list.append(f"-f {self.avoid} ")
            if self.cu and self.factor:
                argument_list.append(f"-fa {self.factor} ")
            if self.template:
                argument_list.append(f"-tem {self.template} ")
            if self.mut:
                argument_list.append(f"-mut {self.mut} ")
            if self.rotamer:
                argument_list.append(f"-rot {self.rotamer} ")
            if self.skip:
                argument_list.append(f"-sk {self.skip} ")
            if self.dihedral_atoms:
                argument_list.append(f"-da {self.dihedral_atoms} ")
                argument_list.append(f"-im {self.improve} ")
            if self.turn:
                argument_list.append(f"-tu {self.turn} ")
            if self.energy_thres:
                argument_list.append(f"-en {self.energy_thres} ")
            if self.conservative:
                argument_list.append(f"-cst {self.conservative} ")
            all_arguments = "".join(argument_list)
            python = f"/gpfs/projects/bsc72/conda_envs/saturated/bin/python -m satumut.simulation {all_arguments}\n"
            lines2.append(python)
            lines.extend(lines2)
            slurm.writelines(lines)

        return self.slurm


def main():
    input_, position, ligchain, ligname, atoms, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, hydrogen, consec, \
    sbatch, steps, dpi, traj, plot_dir, single_mutagenesis, plurizyme_at_and_res, radius, fixed_resids, factor, \
    total_cpus, xtc, cata_dist, template, skip, rotamers, equilibration, log, cpt, improve,  turn, energy_thres, QM, \
    dihedral, box_radius, mut, conservative, profile_with, wild, side_chain_resolution, epochs = parse_args()

    if dir_ and len(input_) > 1:
        dir_ = None
    for inp in input_:
        run = CreateSlurmFiles(inp, ligchain, ligname, atoms, position, cpus, dir_, hydrogen, multiple, pdb_dir, consec,
                               test, cu, seed, nord, steps, dpi, traj, plot_dir, single_mutagenesis,
                               plurizyme_at_and_res, radius, fixed_resids, factor, total_cpus, xtc, cata_dist, template,
                               skip, rotamers, equilibration, log, cpt, improve, turn, energy_thres, QM, dihedral,
                               box_radius, mut, conservative, profile_with, wild, side_chain_resolution, epochs)

        slurm = run.slurm_creation()
        if sbatch:
            call(["sbatch", f"{slurm}"])


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
