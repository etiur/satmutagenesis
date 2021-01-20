# Global imports
import numpy as np
import sys

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

def Neighbourresidues(input_, specific_at_residname, radius=5):
    """
    It gives the list of residues near a specific atom according to a radius
    in a PDB file.

    TODO: It could be implemented with the Biopython library but then the module
    must be loaded

    PARAMETERS
    ----------
    input_ : string
                       PDB file where the specific atoms resides

    specific_at_residname : list of strings
                       PDB atom name of the selected specific atom (spaces in "_"), residue
                       number and residue name

    radius : float
                       Value of the minimum distance between the atom and any of the residues

    Returns
    _______
    Updated_positions : The list of neighbour residues of the specified atom
    """

    Updated_positions = []

    # Get the coordinates of the ester C atom in the ligand to look for the neighbour titrable residues
    PDB_file = open(input_, "rt")
    PDB_file_lines = PDB_file.readlines()
    for line in PDB_file_lines:
        if line[12:16].replace(" ", "_") == specific_at_residname[0] and line[22:26].strip() == specific_at_residname[1] and line[17:20].strip() == specific_at_residname[2]:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])

    # Get the neighbour residues according to the ester carbon in the ligand
    for line in PDB_file_lines:
        if "TER" in line:
            pass
        elif "ATOM" in line:
            x_aux, y_aux, z_aux = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if np.sqrt((x - x_aux) ** 2 + (y - y_aux) ** 2 + (z - z_aux) ** 2) <= float(radius):
                if line[21:22] + ":" + line[22:26].strip() not in Updated_positions:
                    Updated_positions.append(line[21:22] + ":" + line[22:26].strip())
        else:
            pass

    return Updated_positions