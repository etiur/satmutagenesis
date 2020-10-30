from pmx import Model
from pmx.rotamer import mutate
from pmx.rotamer import load_bbdep
import argparse
import os
from helper import map_atom_string

# Argument parsers
def parse_args():
    parser = argparse.ArgumentParser(description="Performs saturated mutagenesis given a PDB file")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, help="Include a chain ID and a position")

    #arguments = vars(parser.parse_args())
    args = parser.parse_args()
    return args.input, args.position

class SaturatedMutagenesis():

    def __init__(self, model, position):
        """
            model (str) path to the PDB file
            position (str) chain ID:position of the residue, for example A:132
        """
        self.model = Model(model)
        self.input = model
        self.coords = position
        self.rotamers = load_bbdep()
        self.residues = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN',
                         'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']
        self.final_pdbs = []
        self.chain = None
        self.chain_id = None
        self.position = None

    def check_coords(self):
        """map the user coordinates with pmx coordinates"""
        if not os.path.exists("pdb_files"):
            os.mkdir("pdb_files")
        self.model.write("pdb_files/original.pdb")
        self.final_pdbs.append("pdb_files/original.pdb")

        after = map_atom_string(self.coords, self.input, "pdb_files/original.pdb")
        self.chain_id = after.split(":")[0]
        self.position = int(after.split(":")[1]) - 1

        for chain_ in self.model.chains:
            if chain_.id == self.chain_id:
                self.chain = chain_

    def generate_pdb(self):
        """
        Generate all the other 19 mutations
        """
        aa_name = self.chain.residues[self.position].resname
        for aa in self.residues:
            if aa != aa_name:
                mutate(self.chain.residues[self.position], aa, self.rotamers)
                output = "{}_{}.pdb".format(aa, self.coords.split(":")[1])
                self.model.write("pdb_files/{}".format(output))
                self.final_pdbs.append("pdb_files/{}".format(output))

        return self.final_pdbs

    def insert_atomtype(self):
        """modifies the pmx PDB files to include the atom type"""
        # read in user input
        with open(self.input, "r") as initial:
            initial_lines = initial.readlines()

        # read in preprocessed input
        for prep_pdb in self.final_pdbs:
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

                elif (line.startswith("HETATM") or line.startswith("ATOM")) and line[
                    21].strip() == self.chain_id.strip() and line[
                                                     22:26].strip() == str(self.position + 1):
                    atom_name = line[12:16].strip()
                    atom_type = "           {}  \n".format(atom_name[0])
                    prep_lines[ind] = line.strip("\n") + atom_type

            with open(prep_pdb, "w") as prep:
                prep.writelines(prep_lines)

def generate_mutations(input, position):
    """
        input (str) - Input pdb to be used to generate the mutations
        position (str) - chain ID:position of teh residue, for example A:139
        start (int) - The position of the first residue in the PDB file
    """
    run = SaturatedMutagenesis(input, position)
    run.check_coords()
    final_pdbs = run.generate_pdb()
    run.insert_atomtype()

    return final_pdbs
    
def main():
    input, position = parse_args()
    output = generate_mutations(input, position)

    return output

if __name__ == "__main__": 
    #Run this if this file is executed from command line but not if is imported as API
    all_pdbs = main()
