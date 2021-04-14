"""
This script is used to generate the yaml files for pele platform
"""

import argparse
import os
from helper import map_atom_string
import glob
from os.path import dirname


def parse_args():
    parser = argparse.ArgumentParser(description="Generate running files for PELE")
    # main required arguments
    parser.add_argument("--folder", required=True,
                        help="An iterable of the path to different pdb files, a name of the folder with the pdbs")
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
    parser.add_argument("-s", "--seed", required=False, default=12345, type=int,
                        help="Include the seed number to make the simulation reproducible")
    parser.add_argument("-st", "--steps", required=False, type=int, default=1000,
                        help="The number of PELE steps")
    parser.add_argument("-x", "--xtc", required=False, action="store_true",
                        help="Change the pdb format to xtc")
    parser.add_argument("-tem", "--template", required=False, nargs="+",
                        help="Path to external forcefield templates")
    parser.add_argument("-rot", "--rotamers", required=False, nargs="+",
                        help="Path to external rotamers templates")
    parser.add_argument("-sk", "--skip", required=False,
                        help="skip the processing of ligands by PlopRotTemp")
    args = parser.parse_args()

    return [args.folder, args.ligchain, args.ligname, args.atoms, args.cpus_per_mutant, args.test, args.polarize_metals,
            args.seed, args.nord, args.steps, args.polarization_factor, args.total_cpus, args.xtc, args.template,
            args.skip, args.rotamers]


class CreateYamlFiles:
    """
    Creates the 2 necessary files for the pele simulations
    """
    def __init__(self, input_path,  ligchain, ligname, atoms, cpus=25,
                 test=False, initial=None, cu=False, seed=12345, nord=False, steps=1000, single=None, factor=None,
                 total_cpus=None, xtc=False, template=None, skip=None, rotamers=None):
        """
        Initialize the CreateLaunchFiles object

        Parameters
        ___________
        input_path: list[str]
            A list of the path to the mutant pdbs
        ligchain: str
            the chain ID where the ligand is located
        ligname: str
            the residue name of the ligand in the PDB
        atoms: list[str]
            list of atom of the residue to follow, in this format --> chain ID:position:atom name
        cpus: int, optional
            How many cpus do you want to use
        test: bool, optional
            Setting the simulation to test mode
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
        """
        self.input = input_path
        self.ligchain = ligchain
        self.ligname = ligname
        self.atoms = atoms[:]
        self.cpus = cpus
        self.test = test
        self.yaml = None
        self.initial = initial
        self.cu = cu
        self.seed = seed
        self.nord = nord
        self.xtc = xtc
        if single and steps == 1000:
            self.steps = 250
        else:
            self.steps = steps
        self.single = single
        self.factor = factor
        if total_cpus:
            self.total_cpu = total_cpus
        else:
            self.total_cpu = len(self.input) * self.cpus + 1
        self.template = template
        self.skip = skip
        self.rotamers = rotamers

    def _match_dist(self):
        """
        match the user coordinates to pmx PDB coordinates
        """
        if self.initial:
            for i in range(len(self.atoms)):
                self.atoms[i] = map_atom_string(self.atoms[i], self.initial, self.input[0])
        else:
            pass

    @staticmethod
    def _search_round():
        """
        Looks at which round of the mutation it is
        """
        count = 1
        round_ = "round_1"
        while os.path.exists(round_):
            count += 1
            round_ = "round_{}".format(count)

        return round_

    def input_creation(self):
        """
        create the .yaml input files for PELE
        """
        self._match_dist()
        if self.single:
            folder = self._search_round()
        else:
            folder = "simulations"
        if not os.path.exists("yaml_files"):
            os.mkdir("yaml_files")
        self.yaml = "yaml_files/simulation.yaml"
        with open(self.yaml, "w") as inp:
            lines = ["system: '{}/*.pdb'\n".format(dirname(self.input[0])), "chain: '{}'\n".format(self.ligchain),
                     "resname: '{}'\n".format(self.ligname), "saturated_mutagenesis: true\n",
                     "seed: {}\n".format(self.seed), "steps: {}\n".format(self.steps),
                     "atom_dist:\n"]
            lines_atoms = ["- '{}'\n".format(atom) for atom in self.atoms]
            lines.extend(lines_atoms)
            if self.xtc:
                lines.append("traj: trajectory.xtc\n")
            if not self.nord:
                lines.append("usesrun: true\n")
            lines.append("working_folder: '{}'\n".format(folder))
            if self.test:
                lines.append("test: true\n")
                self.cpus = 2
                self.total_cpu = len(self.input) * self.cpus + 1
            lines2 = ["cpus: {}\n".format(self.total_cpu), "cpus_per_mutation: {}\n".format(self.cpus),
                      "pele_license: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/license'\n"]
            if not self.nord:
                lines2.append("pele_exec: '/gpfs/projects/bsc72/PELE++/mniv/V1.6.1/bin/PELE-1.6.1_mpi'\n")
            else:
                lines2.append("pele_exec: '/gpfs/projects/bsc72/PELE++/nord/V1.6.1/bin/PELE-1.6.1_mpi'\n")
            if self.cu:
                lines2.append("polarize_metals: true\n")
            if self.cu and self.factor:
                lines2.append("polarization_factor: {}\n".format(self.factor))
            if self.template:
                lines2.append("templates:\n")
                for templates in self.template:
                    lines2.append("- '{}'\n".format(templates))
            if self.rotamers:
                lines2.append("rotamers:\n")
                for rotamers in self.rotamers:
                    lines.append("- '{}'\n".format(rotamers))
            if self.skip:
                lines2.append("skip_ligand_prep:\n- '{}'\n".format(self.skip))
            lines.extend(lines2)
            inp.writelines(lines)

        return self.yaml


def create_20sbatch(pdb_files, ligchain, ligname, atoms, cpus=25, test=False, initial=None,
                    cu=False, seed=12345, nord=False, steps=800, single=None, factor=None,
                    total_cpus=None, xtc=False, template=None, skip=None, rotamers=None):
    """
    creates for each of the mutants the yaml and slurm files

    Parameters
    ___________
    pdb_files: str, list[str]
        the directory to the pdbs or a list of the paths to the mutant pdbs
    ligchain: str
        the chain ID where the ligand is located
    ligname: str
        the residue name of the ligand in the PDB
    atoms: list[str]
        list of atom of the residue to follow, in this format --> chain ID:position:atom name
    cpus: int, optional
        how many cpus do you want to use
    test: bool, optional
        Setting the simulation to test mode
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
    Returns
    _______
    yaml: str
        The input file path for the pele_platform
    """
    if type(pdb_files) == str:
        pdb_list = glob.glob("{}/*.pdb".format(pdb_files))
    else:
        pdb_list = pdb_files
    run = CreateYamlFiles(pdb_list, ligchain, ligname, atoms, cpus, test=test,
                          initial=initial, cu=cu, seed=seed, nord=nord, steps=steps, single=single, factor=factor,
                          total_cpus=total_cpus, xtc=xtc, skip=skip, template=template, rotamers=rotamers)
    yaml = run.input_creation()
    return yaml


def main():
    folder, ligchain, ligname, atoms, cpus, test, cu, seed, nord, steps, factor, total_cpus, xtc, template, \
    skip, rotamers = parse_args()
    yaml_files = create_20sbatch(folder, ligchain, ligname, atoms, cpus=cpus, test=test, cu=cu,
                                 seed=seed, nord=nord, steps=steps, factor=factor, total_cpus=total_cpus, xtc=xtc,
                                 skip=skip, template=template,rotamers=rotamers)

    return yaml_files


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    yaml_list = main()
