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
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, nargs="+",
                        help="Include one or more chain IDs and positions --> ID:position")
    parser.add_argument("--multiple", required=False, action="store_true")
    parser.add_argument("--folder", required=False, default="pdb_files", help="The folder for the pdb_files")
    # arguments = vars(parser.parse_args())
    args = parser.parse_args()
    return args.input, args.position, args.multiple, args.folder


class Mutagenesis:

    def __init__(self, model, position, folder="pdb_files"):
        """
        model (str) path to the PDB file
        position (str) chain ID:position of the residue, for example A:132
        folder (str): The folder where the pdbs are written
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

    def mutate(self, residue, new_aa, bbdep, hydrogens=True):
        if len(new_aa) == 1:
            new_aa = _aacids_dic[new_aa]
        phi = residue.get_phi()
        psi = residue.get_psi()
        rotamers = get_rotamers(bbdep, new_aa, phi, psi, residue=residue, full=True, hydrogens=hydrogens)
        new_r = select_best_rotamer(self.model, rotamers)
        self.model.replace_residue(residue, new_r)

    def check_coords(self, mode=0):
        """
        map the user coordinates with pmx coordinates
        mode (0/1): Acts as a switch, 0 if only 1 mutation per PDB, 1 if 2 mutations per PDB
        """
        if not os.path.exists(self.folder):
            os.mkdir(self.folder)
        if not mode:
            self.model.write("{}/original.pdb".format(self.folder))
            self.final_pdbs.append("{}/original.pdb".format(self.folder))

        after = map_atom_string(self.coords, self.input, "{}/original.pdb".format(self.folder))
        self.chain_id = after.split(":")[0]
        self.position = int(after.split(":")[1]) - 1

        for chain_ in self.model.chains:
            if chain_.id == self.chain_id:
                self.chain = chain_

    def saturated_mutagenesis(self, hydrogens=True, mode=0, name=None):
        """
        Generate all the other 19 mutations
        hydrogens (boolean): Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
        mode (0/1): Acts as a switch, 0 if only 1 mutation per PDB, 1 if 2 mutations per PDB
        name (str): Only used when mode set to 1
        """
        aa_init = self.chain.residues[self.position]
        aa_name = self._invert_aa[aa_init.resname]
        for new_aa in self.residues:
            if new_aa != aa_init.resname:
                self.mutate(aa_init, new_aa, self.rotamers, hydrogens=hydrogens)
                # writing into a pdb
                if not mode:
                    output = "{}{}{}.pdb".format(aa_name, self.position + 1, self._invert_aa[new_aa])
                else:
                    output = "{}_{}{}{}.pdb".format(name, aa_name, self.position + 1, self._invert_aa[new_aa])
                self.model.write("{}/{}".format(self.folder, output))
                self.final_pdbs.append("{}/{}".format(self.folder, output))

        return self.final_pdbs

    def single_mutagenesis(self, new_aa, hydrogens=True):
        """
        Create single mutations
        new_aa: The aa to mutate to in 3 letter code or 1 letter code
        hydrogens (boolean): Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
        """
        aa_init = self.chain.residues[self.position]
        aa_name = self._invert_aa[aa_init.resname]
        self.mutate(aa_init, new_aa, self.rotamers, hydrogens=hydrogens)
        # writing into a pdb
        if len(new_aa) == 1:
            new = new_aa
        elif self._invert_aa.get(new_aa):
            new = self._invert_aa[new_aa]
        else:
            raise Exception("Not recognized aminoacid")

        output = "{}{}{}.pdb".format(aa_name, self.position + 1, new)
        file_ = "{}/{}".format(self.folder, output)
        self.model.write(file_)

        return file_

    def insert_atomtype(self, prep_pdb):
        """
        modifies the pmx PDB files to include the atom type
        prep_pdb (file): PDB files to modify
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
        pros = []
        if file_list:
            self.final_pdbs = file_list
        for prep_pdb in self.final_pdbs:
            p = Process(target=self.insert_atomtype, args=(prep_pdb,))
            p.start()
            pros.append(p)
        for p in pros:
            p.join()


def generate_mutations(input_, position, hydrogens=True, multiple=False, folder="pdb_files"):
    """
    To generate up to 2 mutations per pdb
    input (str): Input pdb to be used to generate the mutations
    position (list): [chain ID:position] of the residue, for example [A:139,..]
    hydrogens (boolean): Leave it true since it removes hydrogens (mostly unnecessary) but creates an error for CYS
    multiple (boolean): Specify if to mutate 2 positions at the same pdb
    """
    count = 0
    pdbs = []
    # Perform single saturated mutations
    for mutation in position:
        run = Mutagenesis(input_, mutation, folder)
        if not count and not os.path.exists("pdb_files/original.pdb"):
            run.check_coords()
        else:
            run.check_coords(mode=1)
        final_pdbs = run.saturated_mutagenesis(hydrogens=hydrogens)
        pdbs.extend(final_pdbs)
        run.accelerated_insert()
        # Mutate in a second position for each of the single mutations
        if multiple and not count and len(position) == 2:
            for files in final_pdbs:
                name = basename(files)
                if name != "original.pdb":
                    name = name.replace(".pdb", "")
                    run_ = Mutagenesis(files, position[1], folder)
                    run_.check_coords(mode=1)
                    final_pdbs_2 = run_.saturated_mutagenesis(hydrogens=hydrogens, mode=1, name=name)
                    pdbs.extend(final_pdbs_2)
                    run_.accelerated_insert()

        count += 1

    return pdbs


def main():
    input_, position, multiple, folder = parse_args()
    output = generate_mutations(input_, position, multiple, folder)

    return output


if __name__ == "__main__":
    # Run this if this file is executed from command line but not if is imported as API
    all_pdbs = main()
