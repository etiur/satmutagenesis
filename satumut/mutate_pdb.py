"""
This script used pmx to mutate the residues within proteins
"""

from pmx import Model
from pmx.rotamer import load_bbdep
import argparse
import os
from helper import map_atom_string
from pmx.library import _aacids_dic
from pmx.rotamer import get_rotamers, select_best_rotamer
from os.path import basename
from multiprocessing import Process


# Argument parsers
def parse_args():
    parser = argparse.ArgumentParser(description="Performs saturated mutagenesis given a PDB file")
    # main required arguments
    parser.add_argument("-i","--input", required=True, help="Include PDB file's path")
    parser.add_argument("-p","--position", required=True, nargs="+",
                        help="Include one or more chain IDs and positions -> Chain ID:position")
    parser.add_argument("-m","--multiple", required=False, action="store_true",
                        help="if you want to mutate 2 residue in the same pdb")
    parser.add_argument("-h","--hydrogen", required=False, action="store_false", help="leave it to default")
    parser.add_argument("-co","--consec", required=False, action="store_true",
                        help="Consecutively mutate the PDB file for several rounds")
    parser.add_argument("-pd","--pdb_dir", required=False, default="pdb_files",
                        help="The name for the mutated pdb folder")
    # arguments = vars(parser.parse_args())
    args = parser.parse_args()
    return args.input, args.position, args.hydrogen, args.multiple, args.pdb_dir, args.consec


class Mutagenesis:
    """
    To perform mutations on PDB files
    """
    def __init__(self, model, position, folder="pdb_files", consec=False):
        """
        Initialize the Mutagenesis object

        Parameters
        ___________
        model: str
           path to the PDB file
        position: str
           chain ID:position of the residue, for example A:132
        folder: str
           The folder where the pdbs are written
        consec: bool
           If this is the second round of mutation
        """
        self.model = Model(model)
        self.input = model
        self.coords = position
        self.rotamers = load_bbdep()
        self.residues = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN',
                         'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']
        self.final_pdbs = []
        self.chain = None
        self.position = None
        self._invert_aa = {v: k for k, v in _aacids_dic.items()}
        self.folder = folder
        self.chain_id = None
        self.consec = consec

    def mutate(self, residue, new_aa, bbdep, hydrogens=True):
        """
        Mutate the wild type residue to a new residue

        Parameters
        ___________
        residue: pmx object
            The residue has to be a pmx object
        new_aa: str
            A 3 letter or 1 letter code to represent the new residue
        bbdep:
            A database that can be interpreted by pmx
        hydrogens: bool, optional
            A boolean, leave it to True because False cause problems with cysteine
        """
        if len(new_aa) == 1:
            new_aa = _aacids_dic[new_aa]
        phi = residue.get_phi()
        psi = residue.get_psi()
        rotamers = get_rotamers(bbdep, new_aa, phi, psi, residue=residue, full=True, hydrogens=hydrogens)
        new_r = select_best_rotamer(self.model, rotamers)
        self.model.replace_residue(residue, new_r)

    def _check_coords(self):
        """
        map the user coordinates with pmx coordinates
        """
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)
        if not os.path.exists("{}/original.pdb".format(self.folder)):
            self.model.write("{}/original.pdb".format(self.folder))
            self.final_pdbs.append("{}/original.pdb".format(self.folder))
        if self.consec:
            self.final_pdbs.remove("{}/original.pdb")

        after = map_atom_string(self.coords, self.input, "{}/original.pdb".format(self.folder))
        self.chain_id = after.split(":")[0]
        self.position = int(after.split(":")[1]) - 1

        for chain_ in self.model.chains:
            if chain_.id == self.chain_id:
                self.chain = chain_

    def saturated_mutagenesis(self, hydrogens=True):
        """
        Generate all the other 19 mutations

        Parameters
        ___________
        hydrogens: bool, optional
            Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
        mode: 0/1, optional
            Acts as a switch, 0 if only 1 mutation per PDB, 1 if 2 mutations per PDB

        Returns
        _______
         final_pdbs: list[path]
            A list of the new files
        """
        self._check_coords()
        aa_init_resname = self.chain.residues[self.position].resname
        aa_name = self._invert_aa[aa_init_resname]
        for new_aa in self.residues:
            if new_aa != aa_init_resname:
                self.mutate(self.chain.residues[self.position], new_aa, self.rotamers, hydrogens=hydrogens)
                # writing into a pdb
                if self.consec:
                    name = basename(self.input).replace("pdb", "")
                    output = "{}_{}{}{}.pdb".format(name, aa_name, self.position + 1, self._invert_aa[new_aa])
                else:
                    output = "{}{}{}.pdb".format(aa_name, self.position + 1, self._invert_aa[new_aa])

                self.model.write("{}/{}".format(self.folder, output))
                self.final_pdbs.append("{}/{}".format(self.folder, output))

        return self.final_pdbs

    def single_mutagenesis(self, new_aa, hydrogens=True):
        """
        Create single mutations

        Parameters
        ___________
        new_aa: str
            The aa to mutate to, in 3 letter code or 1 letter code
        hydrogens: bool, optional
            Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
        mode: 0/1, optional
            0 if it is just 1 mutation per PDB, 1 if there are more than one mutations

        Returns
        ______
        file_: str
            The name of the new pdb file
        """
        self._check_coords()
        aa_init_resname = self.chain.residues[self.position].resname
        aa_name = self._invert_aa[aa_init_resname]
        self.mutate(self.chain.residues[self.position], new_aa, self.rotamers, hydrogens=hydrogens)
        # writing into a pdb
        if len(new_aa) == 1:
            new = new_aa
        elif self._invert_aa.get(new_aa):
            new = self._invert_aa[new_aa]
        else:
            raise Exception("Aminoacid not recognized")

        if self.consec:
            name = basename(self.input).replace("pdb", "")
            output = "{}_{}{}{}.pdb".format(name, aa_name, self.position + 1, new)
        else:
            output = "{}{}{}.pdb".format(aa_name, self.position + 1, new)

        file_ = "{}/{}".format(self.folder, output)
        self.model.write(file_)
        self.insert_atomtype(file_)

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
                coords = line[30:54].split()
                for linex in initial_lines:
                    if linex[30:54].split() == coords:
                        prep_lines[ind] = line.strip("\n") + linex[66:81]
                        break

            elif (line.startswith("HETATM") or line.startswith("ATOM")) and line[
                21].strip() == self.chain_id.strip() and line[
                                                         22:26].strip() == str(self.position + 1):

                atom_name = line[12:16].strip()
                if atom_name[0].isalpha():
                    atom_type = "           {}  \n".format(atom_name[0])
                else:
                    atom_type = "           {}  \n".format(atom_name[1])

                prep_lines[ind] = line.strip("\n") + atom_type

        # rewrittes the files now with the atom type
        with open(prep_pdb, "w") as prep:
            prep.writelines(prep_lines)

    def accelerated_insert(self, file_list=None):
        """
        Paralelizes the insert atomtype function

        Parameters
        ___________
        file_list: list[path]
            optional if you want to include another list
        """
        pros = []
        if file_list:
            self.final_pdbs = file_list
        for prep_pdb in self.final_pdbs:
            p = Process(target=self.insert_atomtype, args=(prep_pdb,))
            p.start()
            pros.append(p)
        for p in pros:
            p.join()


def generate_mutations(input_, position, hydrogens=True, multiple=False, pdb_dir="pdb_files", consec=False):
    """
    To generate up to 2 mutations per pdb

    Parameters
    ___________
    input_: str
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

    Returns
    ________
    pdbs: list[paths]
        The list of all generated pdbs' path
    """
    pdbs = []
    # Perform single saturated mutations
    for mutation in position:
        run = Mutagenesis(input_, mutation, pdb_dir, consec)
        final_pdbs = run.saturated_mutagenesis(hydrogens=hydrogens)
        pdbs.extend(final_pdbs)
        run.accelerated_insert()
        # Mutate in a second position for each of the single mutations
        if multiple and len(position) == 2:
            for files in final_pdbs:
                name = basename(files).replace(".pdb", "")
                if name != "original.pdb":
                    run_ = Mutagenesis(files, position[1], pdb_dir, consec)
                    final_pdbs_2 = run_.saturated_mutagenesis(hydrogens=hydrogens)
                    pdbs.extend(final_pdbs_2)
                    run_.accelerated_insert()

    return pdbs


def main():
    input_, position, hydrogen, multiple, pdb_dir, consec = parse_args()
    output = generate_mutations(input_, position, hydrogen, multiple, pdb_dir, consec)

    return output


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    all_pdbs = main()
