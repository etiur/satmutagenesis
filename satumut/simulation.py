"""
This script is used to create and control the simulations
"""
import argparse
from mutate_pdb import generate_mutations
from pele_files import create_20sbatch
from subprocess import Popen
from os.path import abspath, basename
import os
import logging
import time
from analysis import consecutive_analysis
from helper import Neighbourresidues


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("-i","--input", required=True, help="Include PDB file's path")
    parser.add_argument("-p","--position", required=False, nargs="+",
                        help="Include one or more chain IDs and positions -> Chain ID:position")
    parser.add_argument("-lc","--ligchain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("-ln","--ligname", required=True, help="The ligand residue name")
    parser.add_argument("-at","--atoms", required=True, nargs="+",
                        help="Series of atoms of the residues to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("--cu", required=False, action="store_true", help="used if there are copper in the system")
    parser.add_argument("-t","--test", required=False, action="store_true", help="Used if you want to run a test before")
    parser.add_argument("-n","--nord", required=False, action="store_true",
                        help="used if LSF is the utility managing the jobs")
    parser.add_argument("-m","--multiple", required=False, action="store_true",
                        help="if you want to mutate 2 residue in the same pdb")
    parser.add_argument("-s","--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    parser.add_argument("-d","--dir", required=False,
                        help="The name of the folder for all the simulations")
    parser.add_argument("-pd","--pdb_dir", required=False, default="pdb_files",
                        help="The name for the mutated pdb folder")
    parser.add_argument("-hy","--hydrogen", required=False, action="store_false", help="leave it to default")
    parser.add_argument("-co","--consec", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")
    parser.add_argument("-sb", "--sbatch", required=False, action="store_false",
                        help="True if you want to lanch the simulation right after creating the slurm file")
    parser.add_argument("-st", "--steps", required=False, type=int, default=700,
                        help="The number of PELE steps")
    parser.add_argument("--dpi", required=False, default=800, type=int,
                        help="Set the quality of the plots")
    parser.add_argument("--box", required=False, default=30, type=int,
                        help="Set how many data points are used for the boxplot")
    parser.add_argument("--traj", required=False, default=10, type=int,
                        help="Set how many PDBs are extracted from the trajectories")
    parser.add_argument("--out", required=False, default="summary",
                        help="Name of the summary file created at the end of the analysis")
    parser.add_argument("--plot", required=False,
                        help="Path of the plots folder")
    parser.add_argument("-an", "--analyse", required=False, choices=("energy", "distance", "both"), default="distance",
                        help="The metric to measure the improvement of the system")
    parser.add_argument("--thres", required=False, default=-0.1, type=float,
                        help="The threshold for the improvement which will affect what will be included in the summary")
    parser.add_argument("-sm","--single_mutagenesis",required=False,
                        help="Specifiy the name of the residue that you want the "
                             "original residue to be mutated to. Both 3 letter "
                             "code and 1 letter code can be used.")
    parser.add_argument("-PR","--plurizyme_at_and_res", required=False,
                        help="Specify the chain ID, residue number and the PDB atom name that"
                             "will set the list of the neighbouring residues for the"
                             "next round. Example: chain ID:position:atom name")
    parser.add_argument("-r","--radius", required=False, default=5.0, type=float,
                        help="The radius around the selected atom to search for the other residues")
    parser.add_argument("-f","--fixed_resids",required=False, default=[], nargs='+',
                        help="Specify the list of residues that you don't want"
                             "to have mutated (Must write the list of residue"
                             "numbers)")
    args = parser.parse_args()

    return [args.input, args.position, args.ligchain, args.ligname, args.atoms,args.cpus, args.test,
            args.cu, args.multiple, args.seed, args.dir, args.nord, args.pdb_dir, args.hydrogen, args.consec,
            args.steps, args.dpi, args.box, args.traj, args.out, args.plot, args.analyse, args.thres,
            args.single_mutagenesis, args.plurizyme_at_and_res, args.radius, args.fixed_resids]


class SimulationRunner:
    """
    A class that configures and runs simulations
    """
    def __init__(self, input_, dir_=None):
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
        self.return_code = []

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

        """
        os.chdir("../")
        if not self.dir:
            base = basename(self.input)
            base = base.replace(".pdb", "")
        else:
            base = basename(self.dir)
        hold = "bla"
        folder = []
        for files in pdb_list:
            name = basename(files).replace(".pdb", "")
            if name != "original" and hold != name[:-1]:
                hold = name[:-1]
                folder.append("{}_mutations/{}\n".format(base, hold))
        dirname = "dirnames_{}.txt".format(base)
        with open(dirname, "w") as txt:
            txt.writelines(folder)

        return dirname

    def submit(self, yaml_list):
        """
        Tries to parallelize the call to the pele_platform

        Parameters
        __________
        yaml_list: list[path]
            A list of paths to the yaml files
        """
        platform = "/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8"
        commands = [["{}".format(platform), "-m", "pele_platform.main", "{}".format(yaml)] for yaml in yaml_list]
        self.proc = [Popen(command) for command in commands]
        start = time.time()
        for p in self.proc:
            p.wait()
        end = time.time()

        # creating a log
        logging.basicConfig(filename='simulation_time.log', level=logging.DEBUG, format='%(asctime)s - %(message)s',
                            datefmt='%d-%b-%y %H:%M:%S')
        logging.info("It took {} to run {} simulations".format(end - start, len(yaml_list)))


def saturated_simulation(input_, position, ligchain, ligname, atoms, cpus=25, dir_=None, hydrogen=True,
                         multiple=False, pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345,
                         nord=False, steps=700, dpi=800, box=30, traj=10, output="summary",
                         plot_dir=None, opt="distance", thres=-0.1):
    """
    A function that uses the SimulationRunner class to run saturated mutagenesis simulations

    Parameters
    __________
    input_: str
        The wild type PDB file path
    position: list[str]
        [chain ID:position] of the residue, for example [A:139,..]
    ligchain: str
        the chain ID where the ligand is located
    ligname: str
        the residue name of the ligand in the PDB
    atoms: list[str]
            list of atom of the residue to follow, in this format --> chain ID:position:atom name
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
    """
    simulation = SimulationRunner(input_, dir_)
    input_ = simulation.side_function()
    pdb_names = generate_mutations(input_, position, hydrogens=hydrogen, multiple=multiple, pdb_dir=pdb_dir,
                                   consec=consec)
    yaml_files = create_20sbatch(ligchain, ligname, atoms, cpus=cpus, test=test, initial=input_,
                                 file_=pdb_names, cu=cu, seed=seed, nord=nord, steps=steps)
    simulation.submit(yaml_files)
    dirname = simulation.pele_folders(pdb_names)
    if not test:
        consecutive_analysis(dirname, dpi, box, traj, output, plot_dir, opt, cpus, thres)


def plurizyme_simulations(input_, ligchain, ligname, atoms, single_mutagenesis, plurizyme_at_and_res,
                          radius=5.0, fixed_resids=[], cpus=25, dir_=None, hydrogen=True,
                          pdb_dir="pdb_files", consec=False, test=False, cu=False, seed=12345,
                          nord=False, steps=700):
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
        Set it to true if there are coppers in the system
    seed: int, optional
        A seed number to make the simulations reproducible
    nord: bool, optional
        True if the system is managed by LSF
    steps: int, optional
        The number of PELE steps
    """
    simulation = SimulationRunner(input_, dir_)
    input_ = simulation.side_function()
    # Using the neighbours search to obtain a list of positions to mutate
    position = Neighbourresidues(input_, plurizyme_at_and_res, radius, fixed_resids)
    pdb_names = generate_mutations(input_, position, hydrogen, pdb_dir=pdb_dir, consec=consec, single=single_mutagenesis)
    yaml_files = create_20sbatch(ligchain, ligname, atoms, cpus=cpus, test=test, initial=input_,
                                 file_=pdb_names, cu=cu, seed=seed, nord=nord, steps=steps)
    simulation.submit(yaml_files)
    dirname = simulation.pele_folders(pdb_names)


def main():
    input_, position, ligchain, ligname, atoms, cpus, test, cu, multiple, seed, dir_, nord, pdb_dir, \
    hydrogen, consec, steps, dpi, box, traj, out, plot_dir, analysis, thres, single_mutagenesis, \
    plurizyme_at_and_res, radius, fixed_resids = parse_args()
    if not plurizyme_at_and_res and not single_mutagenesis:
        # if the other 2 flags are not present perform saturated_simulations
        saturated_simulation(input_, position, ligchain, ligname, atoms, cpus, dir_, hydrogen,
                         multiple, pdb_dir, consec, test, cu, seed, nord, steps, dpi, box, traj, out,
                         plot_dir, analysis, thres)
    else:
        # perform plurizymes simulations
        plurizyme_simulations(input_, ligchain, ligname, atoms, single_mutagenesis, plurizyme_at_and_res,
                          radius, fixed_resids, cpus, dir_, hydrogen, pdb_dir, consec, test, cu, seed, nord, steps)


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    main()

