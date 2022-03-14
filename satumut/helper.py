"""
This script contains helper functions
"""
import Bio.PDB
import logging
from os.path import basename, commonprefix
import itertools
from pathlib import Path


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


def match_dist(dihedral_atoms, input_pdb, wild):
    """
    match the user coordinates to pmx PDB coordinates
    """
    wild = Path(wild)
    topology = wild.parent.parent / f"input/{wild.stem}_processed.pdb"
    atom = dihedral_atoms[:]
    for i in range(len(atom)):
        atom[i] = map_atom_string(atom[i], input_pdb, topology)
    return atom


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


def commonlist(folder_list):
    """
    A function that separates the list into sublists grouping the common elements in one place

    Parameters
    ____________
    folder_list: list[str]
        A list of the folder paths

    Returns
    ________
    new_list: list[list]
        A list of a list of mutations in the same position
    """
    new_list = []
    for key, group in itertools.groupby(folder_list, lambda x: basename(x)[:-1]):
        new_list.append(list(group))
    assert len([item for sublist in new_list for item in sublist]) == len(folder_list)
    return new_list


def neighbourresidues(input_, specific_at_res_chainid, radius=5.0, fixed_resids=()):
    """
    It gives the list of residues near a specific atom according to a radius
    in a PDB file.

    PARAMETERS
    ------------
    input_ : string
                    PDB file where the specific atoms resides
    specific_at_res_chainid : strings
                    Chain_ID:position:atom_name, which will be the center around the search for neighbours.
    radius : float
                    Value of the minimum distance between the atom and any of the residues
    fixed_resids: list of integers
                    List of residue numbers of the residues that the user don't want to mutate

    Return
    -------
    updated_positions : The list of neighbour residues of the specified atom

    """
    specific_at_res_chainid = specific_at_res_chainid.split(":")
    updated_positions = []
    parser = Bio.PDB.PDBParser(QUIET=True)

    # Open the PDB file with the Bio module and get the topology of the desired atom to get the coordinates
    structure = parser.get_structure(input_[:-4], input_)
    try:
        target_residue = structure[0][specific_at_res_chainid[0]][int(specific_at_res_chainid[1])]
        target_atom = target_residue[specific_at_res_chainid[2]]
    except ValueError:
        target_residue = structure[0][specific_at_res_chainid[0]][("H_{}".format(specific_at_res_chainid[1]), 1, " ")]
        target_atom = target_residue[specific_at_res_chainid[2]]

    # Get all atoms of the structure and create an instance for a neighbour search around the desired atom
    atoms = Bio.PDB.Selection.unfold_entities(structure[0], 'A')
    ns = Bio.PDB.NeighborSearch(atoms)

    # Get the close residues to the desired atom by a neighbour search
    close_residues = ns.search(target_atom.coord, radius, level='R')

    # Take the output of the neighbour search with biopython and take the positions of the residues that will be mutated
    for close_res in close_residues:
        if not close_res == target_residue:
            if str(close_res.id[1]) not in fixed_resids and close_res.id[0].isspace():
                updated_positions.append(str(close_res.get_parent().id) + ':' + str(close_res.id[1]))

    return updated_positions


class Log:
    """
    A class to keep log of the output from different modules
    """

    def __init__(self, name):
        """
        Initialize the Log class

        Parameters
        __________
        name: str
            The name of the log file
        """
        self._logger = logging.getLogger(__name__)
        self._logger.handlers = []
        self._logger.setLevel(logging.DEBUG)
        self.fh = logging.FileHandler("{}.log".format(name))
        self.fh.setLevel(logging.DEBUG)
        self.ch = logging.StreamHandler()
        self.ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.fh.setFormatter(formatter)
        self.ch.setFormatter(formatter)
        self._logger.addHandler(self.fh)
        self._logger.addHandler(self.ch)

    def debug(self, messages, exc_info=False):
        """
        It pulls a debug message.

        Parameters
        ----------
        messages : str
            The messages to log
        exc_info : bool, optional
            Set to true to include the exception error message            
        """
        self._logger.debug(messages, exc_info=exc_info)

    def info(self, messages, exc_info=False):
        """
        It pulls an info message.

        Parameters
        ----------
        messages : str
            The messages to log
        exc_info : bool, optional
            Set to true to include the exception error message            
        """
        self._logger.info(messages, exc_info=exc_info)

    def warning(self, messages, exc_info=False):
        """
        It pulls a warning message.

        Parameters
        ----------
        messages : str
            The messages to log
        exc_info : bool, optional
            Set to true to include the exception error message            
        """
        self._logger.warning(messages, exc_info=exc_info)

    def error(self, messages, exc_info=False):
        """
        It pulls a error message.

        Parameters
        ----------
        messages : str
            The messages to log
        exc_info : bool, optional
            Set to true to include the exception error message
        """
        self._logger.error(messages, exc_info=exc_info)

    def critical(self, messages, exc_info=False):
        """
        It pulls a critical message.

        Parameters
        ----------
        messages : str
            The messages to log
        exc_info : bool, optional
            Set to true to include the exception error message
        """
        self._logger.critical(messages, exc_info=exc_info)


def find_log(folder_name, wild=None):
    """
    Find the completed log file and gets the path of the different mutations
    Parameters
    ----------
    folder_name: str
        The name of the folder where the simulations are
    wild: str
        The path to the wildtype simulation

    Returns
    -------
    folder: list[str]
        A list of the path to the different mutants
    original: str
        Path to the wild type
    """
    folder = []
    original = wild
    with open(f"{folder_name}/simulations/completed_mutations.log") as log:
        for paths in log:
            dir_ = paths.split()
            if "original" in dir_[1]:
                original = f"{folder_name}/simulations/{dir_[5]}/output/{dir_[1][:-4]}"
            else:
                folder.append(f"{folder_name}/simulations/{dir_[5]}/output/{dir_[1][:-4]}")

    return folder, original


def weighted_median(df, val, weight):
    """
    Calculates the weighted median of a pandas dataframe

    Parameters
    ----------
    df: Dataframe object
        Pandas dataframe
    val: str
        The column name for the values
    weight: str
        The column name of the weights

    Returns
    -------
    df: Dataframe object

    """
    df_sorted = df.sort_values(val)
    cumsum = df_sorted[weight].cumsum()
    cutoff = df_sorted[weight].sum() / 2.
    return df_sorted[cumsum >= cutoff][val].iloc[0]


def check_completed_log(file_name, wild, plot_dir):
    if isiterable(file_name):
        pele_folders = commonlist(file_name)
    elif Path(file_name).exists():
        folder, wild = find_log(file_name, wild)
        pele_folders = commonlist(folder)
    else:
        raise Exception("Pass a file with the path to the different folders")
    if not plot_dir:
        plot_dir = commonprefix(pele_folders[0])
        plot_dir = list(filter(lambda x: "_mut" in x, plot_dir.split("/")))
        plot_dir = plot_dir[0].replace("_mut", "")

    return pele_folders, plot_dir, wild
