import argparse
from mutate_pdb import generate_mutations
from pele_files import create_20sbatch

def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, help="Include a chain ID and a position starting from 0 -> Chain ID:position")
    parser.add_argument("--chain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--resname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format -> chain ID:position:atom name")

    args = parser.parse_args()
    return args.input, args.position, args.chain, args.resname, args.atom1, args.atom2

def main():
    input, position, chain, resname, atom1, atom2 = parse_args()
    pdb_names = generate_mutations(input, position)
    yaml_files, slurm_files = create_20sbatch(chain, resname, atom1, atom2)

if __name__ == "__main__":
    #Run this if this file is executed from command line but not if is imported as API
    main()
