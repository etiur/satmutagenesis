import argparse
from mutate_pdb import generate_mutations, generate_multiple_mutations
from pele_files import create_20sbatch
from subprocess import call
from os.path import abspath, basename
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Generate the mutant PDB and the corresponding running files")
    # main required arguments
    parser.add_argument("--input", required=True, help="Include PDB file's path")
    parser.add_argument("--position", required=True, nargs="+", help="Include a chain ID and a position -> Chain ID:position")
    parser.add_argument("--chain", required=True, help="Include the chain ID of the ligand")
    parser.add_argument("--resname", required=True, help="The ligand residue name")
    parser.add_argument("--atom1", required=True,
                        help="atom of the residue to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--atom2", required=True,
                        help="atom of the ligand to follow in this format -> chain ID:position:atom name")
    parser.add_argument("--cpus", required=False, default=24, type=int,
                       help="Include the number of cpus desired")
    parser.add_argument("--test", required=False, action="store_true")
    parser.add_argument("--cu", required=False, action="store_true")
    parser.add_argument("--multiple", required=False, action="store_true")

    args = parser.parse_args()

    return args.input, args.position, args.chain, args.resname, args.atom1, args.atom2, args.cpus, args.test, args.cu, \
           args.multiple

def submit(slurm_folder):
    """Given a folder submits the job to the supercomputer"""
    for file in slurm_folder:
        call(["sbatch", "{}".format(file)])

def main():
    input_, position, chain, resname, atom1, atom2, cpus, test, cu, multiple = parse_args()
    input_ = abspath(input_)
    base = basename(input_)
    base = base.replace(".pdb", "")
    if not os.path.exists("mutations_{}".format(base)):
        os.mkdir("mutations_{}".format(base))
    os.chdir("mutations_{}".format(base))
    if multiple:
        pdb_names = generate_multiple_mutations(input_, position, hydrogens=True)
    else:
        pdb_names = generate_mutations(input_, position, hydrogens=True)

    slurm_files = create_20sbatch(chain, resname, atom1, atom2, cpus=cpus, test=test, initial=input_, file_list=pdb_names, cu=cu)
    #submit(slurm_files)
    os.chdir("../")

if __name__ == "__main__":
    #Run this if this file is executed from command line but not if is imported as API
    main()
