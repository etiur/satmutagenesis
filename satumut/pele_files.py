"""
This script is used to generate the yaml files for pele platform
"""

import argparse
from .helper import map_atom_string
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Generate running files for PELE")
    # main required arguments
    parser.add_argument("--folder", required=True,
                        help="An iterable of the path to different pdb files, a name of the folder with the pdbs")
    parser.add_argument("-lc", "--ligchain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("-ln", "--ligname", required=True, help="The ligand residue name")
    parser.add_argument("-at", "--atoms", required=False, nargs="+",
                        help="Series of atoms of the residues to follow by PELE during simulation in this format "
                             "-> chain ID:position:atom name")
    parser.add_argument("-cpm", "--cpus_per_mutant", required=False, default=25, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-tcpus", "--total_cpus", required=False, type=int,
                        help="Include the number of cpus desired")
    parser.add_argument("-po", "--polarize_metals", required=False, action="store_true",
                        help="used if there are metals in the system")
    parser.add_argument("-fa", "--polarization_factor", required=False, type=int,
                        help="The number to divide the charges")
    parser.add_argument("-n", "--nord", required=False, action="store_true",
                        help="used if LSF is the utility managing the jobs")
    parser.add_argument("-s", "--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    parser.add_argument("-st", "--steps", required=False, type=int, default=1000,
                        help="The number of PELE steps")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-e", "--equilibration", required=False, action="store_true",
                        help="Set equilibration")
    parser.add_argument("-tem", "--template", required=False, nargs="+",
                        help="Path to external forcefield templates")
    parser.add_argument("-rot", "--rotamers", required=False, nargs="+",
                        help="Path to external rotamers templates")
    parser.add_argument("-sk", "--skip", required=False, nargs="+",
                        help="skip the processing of ligands by PlopRotTemp")
    parser.add_argument("-l", "--log", required=False, action="store_true",
                        help="write logs")
    parser.add_argument("-co", "--consec", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")
    parser.add_argument("-tu", "--turn", required=False, type=int,
                        help="the round of plurizyme generation, not needed for the 1st round")
    parser.add_argument("--QM", required=False,
                        help="The path to the QM charges")
    parser.add_argument("-br","--box_radius", required=False, type=int,
                        help="Radius of the exploration box")
    parser.add_argument("-scr", "--side_chain_resolution", required=False, type=int, default=10,
                        help="Affects the side chain sampling, the smaller the more accurate")
    parser.add_argument("-ep", "--epochs", required=False, type=int, default=1,
                        help="the number of adaptive epochs to run")
    args = parser.parse_args()

    return [args.folder, args.ligchain, args.ligname, args.atoms, args.cpus_per_mutant, args.polarize_metals,
            args.seed, args.nord, args.steps, args.polarization_factor, args.total_cpus, args.xtc, args.template,
            args.skip, args.rotamers, args.equilibration, args.log, args.consec, args.turn, args.QM, args.box_radius,
            args.side_chain_resolution, args.epochs]


class CreateYamlFiles:
    """
    Creates the 2 necessary files for the pele simulations
    """
    def __init__(self, mutant_list,  ligchain, ligname, atoms=None, cpus=25, initial=None, cu=False, seed=12345, nord=False,
                 steps=1000, single=None, factor=None, total_cpus=None, xtc=False, template=None, skip=None,
                 rotamers=None, equilibration=True, log=False, consec=False, turn=None, input_pdb=None, QM=None,
                 box_radius=None, side_chain_resolution=10, epochs=1):
        """
        Initialize the CreateLaunchFiles object

        Parameters
        ___________
        mutant_list: list[str]
            A list of the path to the mutant pdbs
        ligchain: str
            the chain ID where the ligand is located
        ligname: str
            the residue name of the ligand in the PDB
        atoms: list[str]
            list of atom of the residue to follow, in this format --> chain ID:position:atom name
        cpus: int, optional
            How many cpus do you want to use
        initial: file, optional
            The initial PDB file before the modification by pmx
        cu: bool, optional
            Set it to true if there are metals with more than 2 charges (positive or negative) in the system
        seed: int, optional
            A seed number to make the simulations reproducible
        nord: bool, optional
            True if the system is managed by LSF
        steps: int, optional
            The number of PELE steps
        single: str
            Anything that indicates that we are in purizyme mode
        factor: int, optional
            The number to divide the metal charges
        analysis: bool, optional
            True if you want the analysis by pele
        total_cpus: int, optional
            The total number of cpus, it should be a multiple of the number of cpus
        xtc: bool, optional
            Set to True if you want to change the pdb format to xtc
        template: str, optional
            Path to the external forcefield templates
        skip: str, optional
            Skip the processing of ligands by PlopRotTemp
        rotamers: str: optional
            Path to the external rotamers
        equilibration: bool, optional
            True to include equilibration step before PELE
        log: bool, optional
            True to write log files about pele
        round: int, optional
            The round of plurizymer generation
        input_pdb: str, optional
            The pdb file used as input
        QM: str, optional
            Path to the Qm charges
        box_radius: int, optional
            The radius of the exploration box
        side_chain_resolution: int, optional
            The resolution of the side chain sampling, the smaller the better
        """
        self.mutant_list = mutant_list
        self.ligchain = ligchain
        self.ligname = ligname
        if atoms:
            self.atoms = atoms[:]
        else:
            self.atoms = None
        self.cpus = cpus
        self.yaml = None
        self.initial = initial
        self.cu = cu
        self.seed = seed
        self.nord = nord
        self.xtc = xtc
        if single and steps == 1000:
            self.steps = 400
        else:
            self.steps = steps
        self.single = single
        self.factor = factor
        if total_cpus:
            self.total_cpu = total_cpus
        else:
            if self.single:
                self.total_cpu = len(self.mutant_list) * self.cpus + 1
            else:
                self.total_cpu = 4 * self.cpus + 1
        self.template = template
        self.skip = skip
        self.rotamers = rotamers
        self.equilibration = equilibration
        self.log = log
        self.consec = consec
        self.turn = turn
        self.input_pdb = Path(input_pdb)
        self.qm = QM
        self.box = box_radius
        self.resolution = side_chain_resolution
        self.epochs = epochs

    def _match_dist(self):
        """
        match the user coordinates to pmx PDB coordinates
        """
        if self.initial and self.atoms is not None:
            for i in range(len(self.atoms)):
                self.atoms[i] = map_atom_string(self.atoms[i], self.initial, self.mutant_list[0])
        else:
            pass

    def _search_round(self):
        """
        Looks at which round of the mutation it is
        """
        if self.single and not self.turn:
            folder = Path("round_1")
        elif self.single and self.turn:
            folder = Path(f"round_{self.turn}")
        else:
            count = 1
            folder = Path("simulations")
            if self.consec:
                while folder.exists():
                    count += 1
                    folder = Path(f"simulations_{count}")

        return folder

    def yaml_file(self):
        """
        Generating the corresponding yaml files for each of the mutation rounds
        """
        yaml_folder = Path("yaml_files")
        yaml_folder.mkdir(exist_ok=True)
        yaml = yaml_folder/"simulation.yaml"
        if self.consec or self.turn:
            count = 1
            while yaml.exists():
                yaml = yaml_folder/f"simulation_{count}.yaml"
                count += 1
        return yaml

    def yaml_creation(self):
        """
        create the .yaml files for PELE
        """
        self._match_dist()
        folder = self._search_round()
        self.yaml = self.yaml_file()

        with open(self.yaml, "w") as inp:
            lines = [f"system: '{self.mutant_list[0].parent}/*.pdb'\n", f"chain: '{self.ligchain}'\n",
                     f"resname: '{self.ligname}'\n", "saturated_mutagenesis: true\n",
                     f"seed: {self.seed}\n", f"steps: {self.steps}\n", "use_peleffy: true\n"]
            if self.atoms:
                lines.append("atom_dist:\n")
                lines_atoms = [f"- '{atom}'\n" for atom in self.atoms]
                lines.extend(lines_atoms)
            if self.xtc:
                lines.append("traj: trajectory.xtc\n")
            if not self.nord:
                lines.append("usesrun: true\n")
            if self.turn:
                lines.append(f"working_folder: '{folder}/{self.input_pdb.stem}'\n")
            else:
                lines.append(f"working_folder: '{folder}'\n")
            lines2 = [f"cpus: {self.total_cpu}\n", f"cpus_per_mutation: {self.cpus}\n",
                      "pele_license: '/gpfs/projects/bsc72/PELE++/license'\n"]
            if not self.nord:
                lines2.append("pele_exec: '/gpfs/projects/bsc72/PELE++/mniv/V1.7/bin/PELE-1.7_mpi'\n")
            else:
                lines2.append("pele_exec: '/gpfs/projects/bsc72/PELE++/nord4/V1.7.1/bin/PELE-1.7.1_mpi'\n")
            if self.equilibration:
                lines2.append("equilibration: true\n")
            if self.log:
                lines2.append("log: true\n")
            if self.resolution:
                lines.append(f"sidechain_res: {self.resolution}\n")
            if self.cu:
                lines2.append("polarize_metals: true\n")
            if self.qm:
                lines2.append(f"mae_lig: {self.qm}\n")
            if self.cu and self.factor:
                lines2.append(f"polarization_factor: {self.factor}\n")
            if self.template:
                lines2.append("templates:\n")
                for templates in self.template:
                    lines2.append(f" - '{templates}'\n")
            if self.epochs != 1:
                lines2.append(f"iterations: {self.epochs}\n")
            if self.rotamers:
                lines2.append("rotamers:\n")
                for rotamers in self.rotamers:
                    lines2.append(f" - '{rotamers}'\n")
            if self.skip:
                lines2.append("skip_ligand_prep:\n")
                for skip in self.skip:
                    lines2.append(f" - '{skip}'\n")
            if self.box:
                lines2.append(f"box_radius: {self.box}\n")
            lines.extend(lines2)
            inp.writelines(lines)

        return self.yaml


def create_20sbatch(pdb_files, ligchain, ligname, atoms, cpus=25, initial=None, cu=False, seed=12345, nord=False,
                    steps=1000, single=None, factor=None, total_cpus=None, xtc=False, template=None, skip=None,
                    rotamers=None, equilibration=True, log=False, consec=False, turn=None, input_pdb=None, QM=None,
                    box_radius=None, side_chain_resolution=10, epochs=1):
    """
    creates for each of the mutants the yaml and slurm files

    Parameters
    ___________
    pdb_files: str, list[Path]
        the directory to the pdbs or a list of the Path objects to the mutant pdbs
    ligchain: str
        the chain ID where the ligand is located
    ligname: str
        the residue name of the ligand in the PDB
    atoms: list[str]
        list of atom of the residue to follow, in this format --> chain ID:position:atom name
    cpus: int, optional
        how many cpus do you want to use
    initial: file, optional
        The initial PDB file before the modification by pmx if the residue number are changed
    cu: bool, optional
        Set it to true if there are metals in the system in the system
    seed: int, optional
        A seed number to make the simulations reproducible
    nord: bool, optional
        True if the system is managed by LSF
    steps: int, optional
            The number of PELE steps
    single: str, optional
        Anything that indicates that we are in plurizyme mode
    factor: int, optional
        The number to divide the charges of the metals
    analysis: bool, optional
        True if you want the analysis by pele
    total_cpus: int, optional
        The number of total cpus, it should be a multiple of the number of cpus
    xtc: bool, optional
        Set to True if you want to change the pdb format to xtc
    template: str, optional
        Path to the external forcefield templates
    skip: str, optional
        Skip the processing of ligands by PlopRotTemp
    rotamers: str: optional
            Path to the external rotamers
    equilibration: bool, default=False
        True to include equilibration steps before the simulations
    log: bool, optional
        True to write logs abour PELE steps
    consec: bool, optional
        True if it is the second round of mutation
    turn: int, optional
        The round of the plurizymer generation
    input_pdb: str, optional
        The input pdb file
    QM: str, optional
        The path to the QM charges
    box_radius: int, optional
        The radius of the exploration box
    epochs: int, optional
        The number of adaptive epochs to run

    Returns
    _______
    yaml: str
        The input file path for the pele_platform
    """
    if type(pdb_files) == str:
        pdb_list = Path(pdb_files).glob("*.pdb")
    else:
        pdb_list = pdb_files
    run = CreateYamlFiles(pdb_list, ligchain, ligname, atoms, cpus, initial=initial, cu=cu, seed=seed, nord=nord,
                          steps=steps, single=single, factor=factor, total_cpus=total_cpus, xtc=xtc, skip=skip,
                          template=template, rotamers=rotamers, equilibration=equilibration, log=log, consec=consec,
                          turn=turn, input_pdb=input_pdb, QM=QM, box_radius=box_radius,
                          side_chain_resolution=side_chain_resolution, epochs=epochs)
    yaml = run.yaml_creation()
    return yaml


def main():
    folder, ligchain, ligname, atoms, cpus, cu, seed, nord, steps, factor, total_cpus, xtc, template, \
    skip, rotamers, equilibration, log, consec, turn, QM, box_radius, side_chain_resolution, epochs = parse_args()
    yaml_files = create_20sbatch(folder, ligchain, ligname, atoms, cpus=cpus, cu=cu,
                                 seed=seed, nord=nord, steps=steps, factor=factor, total_cpus=total_cpus, xtc=xtc,
                                 skip=skip, template=template, rotamers=rotamers, equilibration=equilibration, log=log,
                                 consec=consec, turn=turn, QM=QM, box_radius=box_radius,
                                 side_chain_resolution=side_chain_resolution, epochs=epochs)

    return yaml_files


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    yaml_list = main()
