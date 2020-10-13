from pmx import *
from pmx.rotamer import mutate
from pmx.rotamer import load_bbdep
import argparse

# Argument parsers
parser = argparse.ArgumentParser(description="Make predictions")
# main required arguments
parser.add_argument("--input", required=True, help="Include PDB file's path")
parser.add_argument("--position", required=True, help="Include a chain ID and a position starting from 0")
arguments = vars(parser.parse_args())


class saturated_mutagenesis():

    def __init__(self, model, chain_id, position, rotamers):
        self.model = model
        self.chain_id = chain_id
        self.position = position
        self.rotamers = rotamers
        self.residues = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER',
                'ARG', 'THR', 'TRP', 'VAL', 'TYR']
        self.chain = None

    def check_chain(self):
        for chain_ in self.model.chains:
            if chain_.id == self.chain_id:
                self.chain = chain_
        try:
            id(self.chain)

        except NameError:
            print "no {} in the Model".format(self.chain_id)
            raise NameError

    def generate_pdb(self):
        for aa in self.residues:
            if aa != self.chain.residues[self.position].resname:
                mutate(self.chain.residues[self.position], aa, self.rotamers)
                self.model.write("{}_{}.pdb".format(aa, self.position))

model = Model(arguments["input"])
chain_id = arguments["position"].split(":")[0]
position = int(arguments["position"].split(":")[1])
rotamers = load_bbdep()

run = saturated_mutagenesis(model, chain_id, position, rotamers)
run.check_chain()
run.generate_pdb()
