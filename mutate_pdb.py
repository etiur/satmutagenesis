from pmx import Model
from pmx.rotamer import mutate
from pmx.rotamer import load_bbdep
import argparse
import os

# Argument parsers
def parse_args():
    parser = argparse.ArgumentParser(description="Performs saturated mutagenesis given a PDB file")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, help="Include a chain ID and a position")
    parser.add_argument("--start", required=False, type=int, default=1,
                        help="The first residue's position in the PDB file")
    #arguments = vars(parser.parse_args())
    args = parser.parse_args()
    return args.input, args.position, args.start

class SaturatedMutagenesis():

    def __init__(self, model, position, start):
        """
            model (str) path to the PDB file
            position (str) chain ID:position of the residue, for example A:132
        """
        self.model = Model(model)
        self.chain_id = position.split(":")[0]
        self.position = int(position.split(":")[1]) - start
        self.rotamers = load_bbdep()
        self.residues = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN',
                         'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']
        self.chain = None

    def check_chain(self):
        """
            check if the chain provided is in fact in the protein
        """
        for chain_ in self.model.chains:
            if chain_.id == self.chain_id:
                self.chain = chain_
        try:
            id(self.chain)
        except NameError:
            raise NameError("no {} in the Model".format(self.chain_id))

    def generate_pdb(self):
        """
        Generate all the other 19 mutations
        """
        if not os.path.exists("pdb_files"):
            os.mkdir("pdb_files")
        final_pdbs = []
        self.model.write("pdb_files/original.pdb")
        final_pdbs.append("original.pdb")

        aa_name = self.chain.residues[self.position].resname
        for aa in self.residues:
            if aa != aa_name:
                mutate(self.chain.residues[self.position], aa, self.rotamers)
                output = "{}_{}.pdb".format(aa, self.position)
                self.model.write("pdb_files/{}".format(output))
                final_pdbs.append(output)

        return final_pdbs

def generate_mutations(input, position, start):
    """
        input (str) - Input pdb to be used to generate the mutations
        position (str) - chain ID:position of teh residue, for example A:139
        start (int) - The position of the first residue in the PDB file
    """
    run = SaturatedMutagenesis(input, position, start)
    run.check_chain()
    final_pdbs = run.generate_pdb()

    return final_pdbs
    
def main():
    input, position, start = parse_args()
    output = generate_mutations(input, position, start)

    return output

if __name__ == "__main__": 
    #Run this if this file is executed from command line but not if is imported as API
    all_pdbs = main()
