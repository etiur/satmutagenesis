from pmx import Model
from pmx.rotamer import mutate
from pmx.rotamer import load_bbdep
import argparse

# Argument parsers

def parse_args():
    parser = argparse.ArgumentParser(description="Make predictions")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, help="Include a chain ID and a position starting from 0")
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
        self.chain_id = position.split(":")[0]
        self.position = int(position.split(":")[1])
        self.rotamers = load_bbdep()
        self.residues = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER',
                'ARG', 'THR', 'TRP', 'VAL', 'TYR']
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
        aa_name = self.chain.residues[self.position].resname
        for aa in self.residues:
            if aa != aa_name:
                mutate(self.chain.residues[self.position], aa, self.rotamers)
                self.model.write("{}_{}.pdb".format(aa, self.position))


def generate_mutations(input, position):
    """
        input (str) - Input pdb to be used to generate the mutations
        position (str) - chain ID:position of teh residue, for example A:139
    """
    run = SaturatedMutagenesis(input, position)
    run.check_chain()
    run.generate_pdb()
    
def main():
    input, position = parse_args()
    generate_mutations(input, position)

if __name__ == "__main__": 
    #Run this if this file is executed from command line but not if is imported as API
    main()
