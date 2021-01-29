# Global imports
import numpy as np
import sys
import Bio.PDB

def map_atom_string(atom_string, initial_pdb, prep_pdb):
    """
    Maps the chain ID and residue number of the original PDb file to the PDB file after pmx

    Parameters
    ___________
    atom_string: str
        The positions to map
    initial_pdb: str
        The original PDB
    prep_pdb: str
        The changed PDB

    Return
    _______
    after: str
        The new atom string or position
    """
    # read in user input
    with open(initial_pdb, "r") as initial:
        initial_lines = initial.readlines()

    # read in preprocessed input
    with open(prep_pdb, "r") as prep:
        prep_lines = prep.readlines()

    # split the atom string or retrieve from the number
    if len(atom_string.split(":")) == 3:
        chain, resnum, atom_name = atom_string.split(":")
        # extract coordinates from user input
        for i in initial_lines:
            if (i.startswith("HETATM") or i.startswith("ATOM")) and i[21].strip() == chain.strip() and i[
                    22:26].strip() == resnum.strip() and i[12:16].strip() == atom_name.strip():

                coords = i[30:54].split()

                # extract coordinates from preprocessed file
                for p in prep_lines:
                    if p[30:54].split() == coords:
                        new_atom_name = p[12:16].strip()
                        new_resnum = p[22:26].strip()
                        new_chain = p[21].strip()
                        after = "{}:{}:{}".format(new_chain, new_resnum, new_atom_name)
                        break
                break

    else:
        chain, resnum = atom_string.split(":")
        for i in initial_lines:
            if (i.startswith("HETATM") or i.startswith("ATOM")) and i[21].strip() == chain.strip() and i[
                                                                    22:26].strip() == resnum.strip():
                coords = i[30:54].split()

                # extract coordinates from preprocessed file
                for p in prep_lines:
                    if p[30:54].split() == coords:
                        new_resnum = p[22:26].strip()
                        new_chain = p[21].strip()
                        after = "{}:{}".format(new_chain, new_resnum)
                        break
                break

    return after

def isiterable(p_object):
    """
    Test if the parameter is an iterable (not a string or a dict) or a file

    Parameters
    ___________
    p_object: object
        Any object

    Returns
    _______
    True: bool
        Returns true if the conditions are met
    """
    if type(p_object) == str or type(p_object) == dict:
        return False
    try:
        iter(p_object)
    except TypeError:
        return False
    return True

def Neighbourresidues(input_, specific_at_res_chainid, radius=5.0, fixed_resids=[]):
    """
    It gives the list of residues near a specific atom according to a radius
    in a PDB file.

    TODO: It could be implemented with the Biopython library but then the module
    must be loaded

    PARAMETERS
    ----------
    input_ : string
                       PDB file where the specific atoms resides

    specific_at_res_chainid : list of strings
                       PDB atom name of the selected specific atom (spaces in "_"), residue
                       number and residue name

    radius : float
                       Value of the minimum distance between the atom and any of the residues

    fixed_resids: list of integers
                        List of residue numbers of the residues that the user don't want to mutate

    Returns
    _______
    Updated_positions : The list of neighbour residues of the specified atom
    """

    # Updated_positions = []
    #
    # # Get the coordinates of the ester C atom in the ligand to look for the neighbour titrable residues
    # PDB_file = open(input_, "rt")
    # PDB_file_lines = PDB_file.readlines()
    # for line in PDB_file_lines:
    #     if line[12:16].replace(" ", "_") == specific_at_residname[0] and line[22:26].strip() == specific_at_residname[1] and line[17:20].strip() == specific_at_residname[2]:
    #         x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
    #
    # # Get the neighbour residues according to the ester carbon in the ligand
    # for line in PDB_file_lines:
    #     if "TER" in line:
    #         pass
    #     elif "ATOM" in line:
    #         x_aux, y_aux, z_aux = float(line[30:38]), float(line[38:46]), float(line[46:54])
    #         if np.sqrt((x - x_aux) ** 2 + (y - y_aux) ** 2 + (z - z_aux) ** 2) <= float(radius):
    #             if line[21:22] + ":" + line[22:26].strip() not in Updated_positions:
    #                 Updated_positions.append(line[21:22] + ":" + line[22:26].strip())
    #     else:
    #         pass
    #
    # return Updated_positions

    Updated_positions = []
    parser = Bio.PDB.PDBParser(QUIET=True)

    # Open the PDB file with the Bio module and get the topology of the desired atom to get the coordinates
    Structure = parser.get_structure(input_[:-4], input_)
    Target_residue = Structure[0][specific_at_res_chainid[0]][int(specific_at_res_chainid[1])]
    Target_atom = Target_residue[specific_at_res_chainid[2]]

    # Get all atoms of the structure and create an instance for a neighbour search around the desired atom
    Atoms = Bio.PDB.Selection.unfold_entities(Structure[0], 'A')
    ns = Bio.PDB.NeighborSearch(Atoms)

    # Get the close residues to the desired atom by a neighbour search
    Close_residues = ns.search(Target_atom.coord, radius, level = 'R')

    # Take the output of the neighbour search with biopython and take the positions of the residues that will be mutated
    for close_res in Close_residues:
        if not close_res == Target_residue:
            if str(close_res.id[1]) not in fixed_resids:
                Updated_positions.append('A'+':'+str(close_res.id[1]))

    return Updated_positions
