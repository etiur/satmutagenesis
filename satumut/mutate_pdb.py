"""
This script used pmx to mutate the residues within proteins
"""

from pmx import Model # change pmx to modeller https://salilab.org/modeller/wiki/Mutate_model
from pmx.rotamer import load_bbdep
import argparse
from .helper import map_atom_string, Log
from pmx.library import _aacids_ext_amber
from pmx.rotamer import get_rotamers, select_best_rotamer
from multiprocessing import Process
from pathlib import Path
import Bio.Align.substitution_matrices as new_mat


# Argument parsers
def parse_args():
    parser = argparse.ArgumentParser(description="Performs saturated mutagenesis given a PDB file")
    # main required arguments
    parser.add_argument("-i", "--input", required=True, nargs="+", help="Include PDB file's path")
    parser.add_argument("-p", "--position", required=True, nargs="+",
                        help="Include one or more chain IDs and positions -> Chain ID:position")
    parser.add_argument("-m", "--multiple", required=False, action="store_true",
                        help="if you want to mutate 2 residue in the same pdb")
    parser.add_argument("-hy", "--hydrogen", required=False, action="store_false", help="leave it to default")
    parser.add_argument("-co", "--consec", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")
    parser.add_argument("-pd", "--pdb_dir", required=False, default="pdb_files",
                        help="The name for the mutated pdb folder")
    parser.add_argument("-sm", "--single_mutagenesis", required=False,
                        help="Specify the name of the residue that you want the "
                             "original residue to be mutated to. Both 3 letter "
                             "code and 1 letter code can be used. You can even specify the protonated states")
    parser.add_argument("-tu", "--turn", required=False, type=int,
                        help="the round of plurizyme generation, not needed for the 1st round")
    parser.add_argument("-mut", "--mutation", required=False, nargs="+",
                        choices=('ALA', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN',
                                 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR'),
                        help="The aminoacid in 3 letter code")
    parser.add_argument("-cst", "--conservative", required=False, choices=(1, 2), default=None, type=int,
                        help="How conservative should the mutations be, choises are 1 and 2")
    parser.add_argument("-w", "--wild", required=False, default=None,
                        help="The path to the folder where the reports from wild type simulation are")
    parser.add_argument("-mp", "--multiple_pdbs", required=False, action="store_true",
                        help="When you are performing single mutations on different pdbs and want to place them in the"
                             "same folder")

    args = parser.parse_args()
    return [args.input, args.position, args.hydrogen, args.multiple, args.pdb_dir, args.consec, args.single_mutagenesis,
            args.turn, args.mutation, args.conservative, args.wild, args.multiple_pdbs]


class Mutagenesis:
    """
    To perform mutations on PDB files
    """
    def __init__(self, model, position, folder="pdb_files", consec=False, single=None, turn=None, mut=None,
                 conservative=None, multiple=False, initial=None, wild_simulation=None, multiple_input=False):
        """
        Initialize the Mutagenesis object

        Parameters
        ___________
        model: str
           path to the PDB file
        position: str
           chain ID:position of the residue, for example A:132
        folder: str, optional
           The folder where the pdbs are written
        consec: bool, optional
           If this is the second round of mutation
        turn: int, optional
            The round of plurizyme generation
        mut: list[str], optional
            A list of specific mutations
        conservative: int, optional
            How conservative should be the mutations according to Blossum62
        multiple: bool, optional
            Same round but double mutations
        initial: str, optional
            The initial input pdb, used if multiple true to check the coordinates
        wild_simulation: str, optional
            The path to the wild type simulation
        single: str, optional
            The new residue to be mutated to, in 3-letter or 1-letter code
        multiple_input: bool, optional
            If the goal is to mutate several pdbs for the single mutagenesis and want to place them in a single folder
        """
        self.model = Model(str(model))
        self.input = Path(model)
        self.wild = wild_simulation
        if not initial:
            self.initial = self.input
        else:
            self.initial = Path(initial)
        self.multiple_pdbs = multiple_input
        self.coords = position
        self.rotamers = load_bbdep()
        self.final_pdbs = []
        self.position = None
        self._invert_aa = {v: k for k, v in _aacids_ext_amber.items()}
        self._invert_aa["HIS"] = "H"
        self.chain_id = None
        self.folder = Path(folder)
        self.consec = consec
        self.multiple = multiple
        self.log = Log("mutate_errors")
        self.single = single
        self.turn = turn
        self.residues = self.__check_all_and_return_residues(mut, conservative)
        self.initial_connect = self.get_coordinates_from_connect_lines()

    def __check_all_and_return_residues(self, mut, conservative):
        self._check_folders_single_mutations()
        self._check_folders_multiple_mutations()
        self.folder.mkdir(parents=True, exist_ok=True)
        self._check_coords()
        self.aa_init_resname = self.model.residues[self.position].resname
        if not mut and not conservative:
            residues = ['ALA', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN',
                        'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']
        elif mut and not conservative:
            residues = mut
        elif conservative and not mut:
            residues = self._mutation_library_new(conservative)
        return residues

    def _mutation_library_new(self, library=1):
        """
        Determines how conservative should be the mutations in the new versions of biopython

        Parameters
        ___________
        library: int
            Choose between 1 and 2 to configure how conservative should be the mutations
        """
        aa = self._invert_aa[self.aa_init_resname]
        matrix = new_mat.load("BLOSUM62")
        matrix = {k: v for k, v in matrix.items() if "X" not in k and "B" not in k and "Z" not in k and "*" not in k}
        blosum = [key for key in matrix.keys() if aa == key[1] and key.count(aa) < 2 and "P" not in key]
        value = [matrix[x] for x in blosum]
        new_dict = dict(zip([_aacids_ext_amber[x[1]] if x[0] == aa else _aacids_ext_amber[x[0]] for x in blosum], value))
        if library == 1:
            reduced_dict = {k: v for k, v in new_dict.items() if v >= 0}
        elif library == 2:
            reduced_dict = {k: v for k, v in new_dict.items() if v >= -1}

        return reduced_dict.keys()

    def _check_folders_single_mutations(self):
        """
        Check the presence of different folders in single mutagenesis
        """
        if self.turn and self.single and not self.multiple_pdbs:
            if not Path(f"{self.folder.parent/self.folder.name}_1_round_{self.turn}").exists():
                self.folder = Path(f"{self.folder.parent/self.folder.name}_1_round_{self.turn}")
            else:
                files = [x for x in self.folder.parent.glob("*") if self.folder.name in x.name]
                files.sort(key=lambda x: int(x.name.replace(f"_round_{self.turn}", "").split("_")[-1]))
                num = int(files[-1].name.replace(f"_round_{self.turn}", "").split("_")[-1])
                self.folder = Path(f"{self.folder.name}_{num+1}_round_{self.turn}")

    def _check_folders_multiple_mutations(self):
        """
        Check the presence of different folders in saturated mutagenesis
        """
        if self.consec and not self.multiple:
            count = 1
            self.folder = Path("next_round_1")
            while self.folder.exists():
                count += 1
                self.folder = Path(f"next_round_{count}")
        elif self.consec and self.multiple:
            files = [x for x in Path.cwd().glob("*") if "next_round" in x.name or self.folder.name in x.name]
            files.sort(key=lambda x: int(x.name.split("_")[-1]) if x.name.split("_")[-1].isdigit() else -99999)
            self.folder = files[-1]

    def _check_coords(self):
        """
        map the user coordinates with pmx coordinates
        """

        original = self.folder.joinpath("original.pdb")
        if not original.exists():
            self.model.write(str(original))
            self.final_pdbs.append(original)
        after = map_atom_string(self.coords, self.initial, original)
        self.chain_id = after.split(":")[0]
        self.position = int(after.split(":")[1]) - 1
        if self.wild or (self.multiple_pdbs and self.single):
            self.final_pdbs = [x for x in self.final_pdbs if "original" not in x.name]
            original.unlink(missing_ok=True)

    def mutate(self, residue, new_aa, bbdep, hydrogens=True):
        """
        Mutate the wild type residue to a new residue

        Parameters
        ___________
        residue: pmx object
            The residue has to be a pmx object
        new_aa: str
            A 3 letter or 1-letter code to represent the new residue
        bbdep:
            A database that can be interpreted by pmx
        hydrogens: bool, optional
            A boolean, leave it to True because False cause problems with cysteine
        """
        if len(new_aa) == 1:
            new_aa = _aacids_ext_amber[new_aa]
        phi = residue.get_phi(degree=True)
        psi = residue.get_psi(degree=True)
        rotamers = get_rotamers(bbdep, new_aa, phi, psi, residue=residue, full=True, hydrogens=hydrogens)
        new_r = select_best_rotamer(self.model, rotamers)
        self.model.replace_residue(residue, new_r)

    def saturated_mutagenesis(self, hydrogens=True, count=0):
        """
        Generate all the other 19 mutations

        Parameters
        ___________
        hydrogens: bool, optional
            Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS

        Returns
        _______
         final_pdbs: list[path]
            A list of the new files
        """
        aa_name = self._invert_aa[self.aa_init_resname]
        for new_aa in self.residues:
            if new_aa != self.aa_init_resname:
                try:
                    self.mutate(self.model.residues[self.position], new_aa, self.rotamers, hydrogens=hydrogens)
                except KeyError:
                    self.log.error(f"position {self.chain_id}:{self.position+1} has no rotamer in the library so it was skipped", exc_info=True)
                # writing into a pdb
                if self.consec or count == 1:
                    output = Path(f"{self.input.stem}_{aa_name}{self.position+1}{self._invert_aa[new_aa]}.pdb")
                else:
                    output = Path(f"{aa_name}{self.position+1}{self._invert_aa[new_aa]}.pdb")

                self.model.write(str(self.folder/output))
                self.final_pdbs.append(self.folder/output)

        return self.final_pdbs

    def single_mutagenesis(self, hydrogens=True):
        """
        Create single mutations

        Parameters
        ___________
        hydrogens: bool, optional
            Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS

        Returns
        ______
        file_: str
            The name of the new pdb file
        """
        aa_name = self._invert_aa[self.aa_init_resname]
        try:
            self.mutate(self.model.residues[self.position], self.single, self.rotamers, hydrogens=hydrogens)
        except KeyError:
            self.log.error(f"position {self.chain_id}:{self.position+1} has no rotamer in the library so it was skipped",
                           exc_info=True)
        # writing into a pdb
        if len(self.single) == 1:
            new = self.single
        else:
            new = self._invert_aa[self.single]
        if self.turn or self.multiple_pdbs:
            output = Path(f"{self.input.stem}_{aa_name}{self.position+1}{new}.pdb")
        else:
            output = Path(f"{aa_name}{self.position+1}{new}.pdb")

        file_ = self.folder/output
        self.model.write(str(file_))
        self.insert_atomtype(file_)
        self.insert_conect_lines(file_)

        return file_

    def insert_atomtype(self, prep_pdb):
        """
        modifies the pmx PDB files to include the atom type

        Parameters
        ___________
        prep_pdb: path
            PDB files to modify
        """
        # read in user input
        with open(self.input, "r") as initial:
            initial_lines = initial.readlines()

        # read in preprocessed input
        with open(prep_pdb, "r") as prep:
            prep_lines = prep.readlines()

        for ind, line in enumerate(prep_lines):
            if (line.startswith("HETATM") or line.startswith("ATOM")) and (
                    line[21].strip() != self.chain_id.strip() or line[
                                                                 22:26].strip() != str(self.position + 1)):
                coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                for linex in initial_lines:
                    if linex.startswith("HETATM") or linex.startswith("ATOM"):
                        if [float(linex[30:38]), float(linex[38:46]), float(linex[46:54])] == coords:
                            prep_lines[ind] = line.strip("\n") + linex[66:]
                            break

            elif (line.startswith("HETATM") or line.startswith("ATOM")) and line[
                21].strip() == self.chain_id.strip() and line[
                                                         22:26].strip() == str(self.position + 1):

                atom_name = line[12:16].strip()
                if atom_name[0].isalpha():
                    atom_type = f"           {atom_name[0]}  \n"
                else:
                    atom_type = f"           {atom_name[1]}  \n"

                prep_lines[ind] = line.strip("\n") + atom_type

        # rewrittes the files now with the atom type
        with open(prep_pdb, "w") as prep:
            prep.writelines(prep_lines)

    def insert_conect_lines(self, prep_pdb):
        """
        Insert the conect lines to a pdb
        """
        if self.initial_connect:
            # I get both the atom coord and indices from the pmx mutants
            pmx_file_indices = self.get_atom_indices(prep_pdb)
            pmx_all_connect = []
            for x in self.initial_connect:
                # Using the coords of the input pdb I get the atom indices of the pmx file
                # the conect lines are multiple of 5, so if the length is < 5 add white spaces
                pmx_connect = [" "*(5-len(str(pmx_file_indices[y]))) + str(pmx_file_indices[y]) for y in x]
                pmx_connect = f"CONECT{''.join(pmx_connect)}\n"
                pmx_all_connect.append(pmx_connect)
            # read in preprocessed input and write the new connect lines
            with open(prep_pdb, "r") as prep:
                prep_lines = prep.readlines()

            with open(prep_pdb, "w") as prep:
                prep_lines[len(prep_lines)-1:len(prep_lines)-1] = pmx_all_connect
                prep.writelines(prep_lines)

    def accelerated_connect(self):
        """
        Parallelize the insert_conect_lines function
        """
        pros = []
        if self.initial_connect:
            for prep_pdb in self.final_pdbs:
                p = Process(target=self.insert_conect_lines, args=(prep_pdb,))
                p.start()
                pros.append(p)
            for p in pros:
                p.join()

    def accelerated_insert(self):
        """
        Parallelize the insert atomtype function
        """
        pros = []
        for prep_pdb in self.final_pdbs:
            p = Process(target=self.insert_atomtype, args=(prep_pdb,))
            p.start()
            pros.append(p)
        for p in pros:
            p.join()

    def get_coordinates_from_connect_lines(self):
        """
        I transform the atom indices of the input pdb into a list of a list of coordinates of those in CONECT lines
        """
        coords = []
        initial_atom_indices = self.get_atom_indices(self.initial)
        with open(self.initial, "r") as initial:
            initial_lines = initial.readlines()
        connect = [x for x in initial_lines if x.startswith("CONECT")]
        if connect:
            for x in connect:
                x = x.replace("CONECT", "")
                x = x.strip("\n")
                num = len(x) / 5
                if num.is_integer():
                    new_x = [int(x[i * 5:(i * 5) + 5]) for i in range(int(num))] # I get the atom indices of the CONECT
                    new_coords = [initial_atom_indices[y] for y in new_x] # I get the coords using atom indices
                    coords.append(new_coords)

            return coords
        else:
            return []

    def get_atom_indices(self, file):
        """
        I get the indices and coords of the pdb file
        """
        atom_indices = {}
        with open(file, "r") as initial:
            initial_lines = initial.readlines()

        for ind, line in enumerate(initial_lines):
            if line.startswith("HETATM") or line.startswith("ATOM"):
                index = int(line[6:11])
                coords = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atom_indices[coords] = index
                atom_indices[index] = coords

        return atom_indices


def generate_single_mutations(input_, position, single, hydrogens=True, pdb_dir="pdb_files", turn=None, wild=None,
                              multiple_inputs=False):
    """
        To generate different single mutations on a single pdb or the same single mutation to different pdbs

        Parameters
        ___________
        input_: list[str]
            Input pdbs to be used to generate the mutations
        position: list[str] -> if multiple_input provided then each position should match the input
            [chain ID:position] of the residue, for example [A:139,..]
        hydrogens: bool, optional
            Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
        pdb_dir: str, optional
            The name of the folder where the mutated PDB files will be stored
        single: str
            The new residue to mutate the positions to, in 3 letter or 1-letter code
        turn: int, optional
            The round of plurizymer generation
        consec: bool, optional
            If to keep the name of the input pdb
        multiple_inputs: bool, optional
            If the goal is to mutate several pdbs for the single mutagenesis and want to place them in a single folder

        Returns
        ________
        pdbs: list[paths]
            The list of all generated pdbs' path
        """
    pdbs = []
    if not multiple_inputs:
        for mutation in position:
            # If the single_mutagenesis flag is used, execute this
            single = single.upper()
            run = Mutagenesis(input_[0], mutation, pdb_dir, single=single, turn=turn, wild_simulation=wild,
                              multiple_input=multiple_inputs)
            mutant = run.single_mutagenesis(hydrogens)
            pdbs.append(mutant)

        if not wild:
            ori = run.folder/"original.pdb"
            run.insert_atomtype(ori)
            run.insert_conect_lines(ori)
            pdbs.append(run.folder/"original.pdb")
    else:
        for num, pdb in enumerate(input_):
            single = single.upper()
            run = Mutagenesis(pdb, position[num], pdb_dir, single=single, turn=turn, wild_simulation=wild,
                              multiple_input=multiple_inputs)
            mutant = run.single_mutagenesis(hydrogens)
            pdbs.append(mutant)

    return pdbs


def generate_saturated_mutations(input_, position, hydrogens=True, multiple=False, pdb_dir="pdb_files", consec=False,
                                 mut=None, conservative=None, wild=None):
    """
    To generate up to 2 mutations per pdb

    Parameters
    ___________
    input_: list[str]
        Input pdb to be used to generate the mutations
    position: list[str]
        [chain ID:position] of the residue, for example [A:139,..]
    hydrogens: bool, optional
        Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
    multiple: bool, optional
        Specify if to mutate 2 positions at the same pdb
    pdb_dir: str, optional
        The name of the folder where the mutated PDB files will be stored
    consec: bool, optional
        Consecutively mutate the PDB file for several rounds
    mut: list[str]
        A list of mutations to perform
    conservative: int, optional
        How conservative should be the mutations according to Blossum62

    Returns
    ________
    pdbs: list[paths]
        The list of all generated pdbs' path
    """
    pdbs = []
    # Perform single saturated mutations
    count = 0
    for mutation in position:
        if multiple and count == 1:
            run = Mutagenesis(input_[0], mutation, pdb_dir, consec, mut=mut, conservative=conservative, multiple=multiple,
                              wild_simulation=wild)
        else:
            run = Mutagenesis(input_[0], mutation, pdb_dir, consec, mut=mut, conservative=conservative, wild_simulation=wild)
        # run saturated mutagenesis
        final_pdbs = run.saturated_mutagenesis(hydrogens=hydrogens)
        pdbs.extend(final_pdbs)
        run.accelerated_insert()
        run.accelerated_connect()
        count += 1
        # Mutate in a second position for each of the 20 single mutations
        if multiple and count == 1:
            for files in final_pdbs:
                if files.stem != "original":
                    run_ = Mutagenesis(files, position[1], pdb_dir, consec, conservative=conservative, mut=mut,
                                       multiple=multiple, initial=input_)
                    final_pdbs_2 = run_.saturated_mutagenesis(hydrogens=hydrogens, count=count)
                    pdbs.extend(final_pdbs_2)
                    run_.accelerated_insert()
                    run_.accelerated_connect()

    return pdbs


def generate_mutations(input_, position, hydrogen=True, multiple=False, pdb_dir="pdb_files", consec=False, single=None,
                       turn=None, mut=None, conservative=None, wild=None, multiple_inputs=False):
    """
    A function that combines both previous functions
    """
    if single:
        pdbs = generate_single_mutations(input_, position, single, hydrogen, pdb_dir, turn, wild, multiple_inputs)
    else:
        pdbs = generate_saturated_mutations(input_, position, hydrogen, multiple, pdb_dir, consec, mut, conservative,
                                            wild)
    return pdbs


def main():
    input_, position, hydrogen, multiple, pdb_dir, consec, single_mutagenesis, turn, mut, conservative, wild,\
    multiple_inputs= parse_args()
    output = generate_mutations(input_, position, hydrogen, multiple, pdb_dir, consec, single_mutagenesis, turn, mut,
                                conservative, wild, multiple_inputs)

    return output


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    all_pdbs = main()
