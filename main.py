import argparse
from mutate_pdb import generate_mutations
from pele_files import create_20sbatch
from subprocess import call
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, help="Include a chain ID and a position -> Chain ID:position")
    parser.add_argument("--chain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--resname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=24, type=int,
                       help="Include the number of cpus desired")
    parser.add_argument("--test", required=False, action="store_true")

    args = parser.parse_args()

    return args.input, args.position, args.chain, args.resname, args.atom1, args.atom2, args.cpus, args.test

def submit(slurm_folder):
    """Given a folder submits the job to the supercomputer"""
    for file in glob.glob("{}/*".format(slurm_folder)):
        call(["sbatch", "{}".format(file)])

def main():
    input, position, chain, resname, atom1, atom2, cpus, test = parse_args()
    pdb_names = generate_mutations(input, position, hydrogens=True)
    yaml_files, slurm_files = create_20sbatch(chain, resname, atom1, atom2, cpus=cpus, test=test, initial=input)
    submit("slurm_files")

if __name__ == "__main__":
    #Run this if this file is executed from command line but not if is imported as API
    main()
