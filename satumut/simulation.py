"""
This script is used to create and control the simulations
"""
import argparse
from mutate_pdb import generate_mutations
from pele_files import create_20sbatch
from subprocess import call, Popen
from os.path import abspath, basename
import os
import time
from analysis import consecutive_analysis
from helper import neighbourresidues, Log


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("-i", "--input", required=True, help="Include PDB file's path")
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
                             "to have mutated (Must write the list of residue"
                             "numbers)")
    parser.add_argument("-re", "--restart", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-cd", "--catalytic_distance", required=False, default=3.5, type=float,
                        help="The distance considered to be catalytic")
    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atoms, args.cpus_per_mutant, args.test,
            args.polarize_metals, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen,
            args.consec, args.steps, args.dpi, args.box, args.trajectory, args.out, args.plot, args.analyse, args.thres,
            args.single_mutagenesis, args.plurizyme_at_and_res, args.radius, args.fixed_resids,
            args.polarization_factor, args.total_cpus, args.restart, args.xtc, args.catalytic_distance]


class SimulationRunner:
    """
    A class that configures and runs simulations
    """

    def __init__(self, input_, dir_=None, single=None, yaml=None):
        """
        Initialize the Simulation Runner class

        Parameters
        ___________
        input_: str
            The path to the PDB file
        dir_: str, optional
            The name of the directory for the simulations to run and the outputs to be stored
        """

        self.input = input_
        self.proc = None
        self.dir = dir_
        self.single = single
        self.log = Log("simulation_time")
        self.yaml = yaml

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
            base = self.input.replace(".pdb", "")
        else:
            base = self.dir.replace(".pdb", "")
        if not os.path.exists("{}_mut".format(base)):
            os.makedirs("{}_mut".format(base))
        os.chdir("{}_mut".format(base))

        return self.input

    def pele_folders(self):
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
        folder = []
        if not self.single:
            with open("{}_mutations/simulations/completed_mutations.log".format(base)) as log:
                for paths in log:
                    dir_ = paths.split()
                    if "original" in dir_[1]:
                        original = "{}_mutations/simulations/{}/output/{}".format(base, dir_[5], dir_[1][:-4])
                    else:
                        folder.append("{}_mutations/simulations/{}/output/{}".format(base, dir_[5], dir_[1][:-4]))
            return folder, original

    def submit(self, yaml):
        """
        Tries to parallelize the call to the pele_platform

        Parameters
        __________
        yaml_list: list[path]
            A list of paths to the yaml files
        """
        command = ["python", "-m", "pele_platform.main", "{}".format(yaml)]
        start = time.time()
        retun_code = call(command, close_fds=False)
        end = time.time()
        # creating a log
        self.log.info("It took {} to run the simulation with return code {}".format(end - start, retun_code))


def saturated_simulation(input_, ligchain, ligname, atoms, position=None, cpus=25, dir_=None, hydrogen=True,
                         multiple=False, pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345,
                         nord=False, steps=800, dpi=800, box=30, traj=10, output="summary",
                         plot_dir=None, opt="distance", thres=-0.1, factor=None, plurizyme_at_and_res=None,
                         radius=5.0, fixed_resids=(), total_cpus=None, restart=False, cata_dist=3.5, xtc=False):
    """
    A function that uses the SimulationRunner class to run saturated mutagenesis simulations

    Parameters
    __________
    input_: str
        The wild type PDB file path
    ligchain: str
        the chain ID where the ligand is located
    ligname: str
        the residue name of the ligand in the PDB
    atoms: list[str]
            list of atom of the residue to follow, in this format --> chain ID:position:atom name
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
    factor: int, optional
        The number to divide the metal charges
    plurizyme_at_and_res: str
        Chain_ID:position:atom_name, which will be the center around the search for neighbours
    radius: float, optional
        The radius for the neighbours search
    fixed_resids. list[position_num]
        A list of residues positions to avoid mutating
    total_cpus: int, optional
        Set the total number of cpus, it should be a multiple of the number of cpus
    restart: bool, optional
        True if the simulation has already run once
    cata_dist: float, optional
        The catalytic distance
    xtc: bool, optional
        Set to True if you want to change the pdb format to xtc
    """
    simulation = SimulationRunner(input_, dir_)
    input_ = simulation.side_function()
    if not position and plurizyme_at_and_res:
        position = neighbourresidues(input_, plurizyme_at_and_res, radius, fixed_resids)
    if not restart:
        pdb_names = generate_mutations(input_, position, hydrogens=hydrogen, multiple=multiple, pdb_dir=pdb_dir,
                                   consec=consec)
        yaml = create_20sbatch(pdb_names, ligchain, ligname, atoms, cpus=cpus, test=test, initial=input_,
                               cu=cu, seed=seed, nord=nord, steps=steps, factor=factor,
                               total_cpus=total_cpus, xtc=xtc)
    else:
        yaml = "yaml_files/simulation.yaml"
        with open(yaml, "r") as yml:
            if "adaptive_restart: true\n" not in  yml.readlines():
                with open(yaml, "a") as yam:
                    yam.write("adaptive_restart: true\n")

    simulation.submit(yaml)
    dirname, original = simulation.pele_folders()
    if not test:
        if dir_ and not plot_dir:
            plot_dir = dir_
        consecutive_analysis(dirname, original, dpi, box, traj, output, plot_dir, opt, cpus, thres, cata_dist, xtc)


def plurizyme_simulation(input_, ligchain, ligname, atoms, single_mutagenesis, plurizyme_at_and_res,
                         radius=5.0, fixed_resids=(), cpus=30, dir_=None, hydrogen=True,
                         pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345,
                         nord=False, steps=250, factor=None, total_cpus=None, xtc=False):
    """
    Run the simulations for the plurizyme's projct which is based on single mutations

    Parameters
    __________
    input_: str
        The wild type PDB file path
    ligchain: str
        the chain ID where the ligand is located
    ligname: str
        the residue name of the ligand in the PDB
    atoms: list[str]
        list of atom of the residue to follow, in this format --> chain ID:position:atom name
    single_mutagenesis: str
        The new residue to mutate the positions to, in 3 letter or 1 letter code
    plurizyme_at_and_res: str
        Chain_ID:position:atom_name, which will be the center around the search for neighbours
    radius: float, optional
        The radius for the neighbours search
    fixed_resids. list[position_num]
        A list of residues positions to avoid mutating
    cpus: int, optional
        how many cpus do you want to use
    dir_: str, optional
        Name of the folder ofr the simulations
    hydrogens: bool, optional
        Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
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
    factor: int, optional
        The number to divide the metal charges
    total_cpus: int, optional
        Set the total number of cpus, it should be a multiple of the number of cpus
    xtc: bool, optional
        Set to True if you want to change the pdb format to xtc
    """
    simulation = SimulationRunner(input_, dir_, single_mutagenesis)
    input_ = simulation.side_function()
    # Using the neighbours search to obtain a list of positions to mutate
    position = neighbourresidues(input_, plurizyme_at_and_res, radius, fixed_resids)
    pdb_names = generate_mutations(input_, position, hydrogen, pdb_dir=pdb_dir, consec=consec,
                                   single=single_mutagenesis)
    yaml = create_20sbatch(pdb_names, ligchain, ligname, atoms, cpus=cpus, test=test, initial=input_,
                                 cu=cu, seed=seed, nord=nord, steps=steps, single=single_mutagenesis,
                                 factor=factor, total_cpus=total_cpus, xtc=xtc)
    simulation.submit(yaml)


def main():
    input_, position, ligchain, ligname, atoms, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, \
    hydrogen, consec, steps, dpi, box, traj, out, plot_dir, analyze, thres, single_mutagenesis, \
    plurizyme_at_and_res, radius, fixed_resids, factor, total_cpus, restart, xtc, cata_dist = parse_args()

    if plurizyme_at_and_res and single_mutagenesis:
        # if the other 2 flags are present perform plurizyme simulations
        plurizyme_simulation(input_, ligchain, ligname, atoms, single_mutagenesis, plurizyme_at_and_res,
                             radius, fixed_resids, cpus, dir_, hydrogen, pdb_dir, consec, test, cu, seed, nord, steps,
                             factor, total_cpus, xtc)
    else:
        # Else, perform saturated mutagenesis
        saturated_simulation(input_, ligchain, ligname, atoms, position, cpus, dir_, hydrogen,
                             multiple, pdb_dir, consec, test, cu, seed, nord, steps, dpi, box, traj, out,
                             plot_dir, analyze, thres, factor, total_cpus, cata_dist, xtc)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()
