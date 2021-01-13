The manual for the Saturated Mutageneis
*****************************************

| The main usage of ``saturated_mutagenesis`` is to mutate a given position within a protein to all the other 19 aminoacids, termed **saturated mutagenesis** and to facilitate their analysis through PELE simulations. 
| As a results it outputs 19 PDBs + 1 PDB for the wildtype and the correspoding files for the PELE simulations in marenostrum or NORD 

Simple explanation
===================
| After the download from the `repository <https://github.com/etiur/saturated_mutagenesis>`_ you can readily use the main.py through the command line to generate the different mutations and simulation files.
| Let's see the necessary arguments

.. code-block:: bash

    $ python saturated_mutagenesis.main --help

.. code-block::

    usage: main.py [-h] --input INPUT --position POSITION [POSITION ...]
               --ligchain LIGCHAIN --ligname LIGNAME --atom1 ATOM1 --atom2
               ATOM2 [--cpus CPUS] [--cu] [--test] [--nord] [--multiple]
               [--seed SEED] [--dir DIR] [--pdb_dir PDB_DIR] [--hydrogen]
               [--consec]

    Generate the mutant PDB and the corresponding running files

    optional arguments:
        -h, --help            show this help message and exit
        --input INPUT         Include PDB file's path
        --position POSITION [POSITION ...]
                        Include one or more chain IDs and positions -> chain ID:position
        --ligchain LIGCHAIN   Include the chain ID of the ligand
        --ligname LIGNAME     The ligand residue name
        --atom1 ATOM1         atom of the residue to follow in this format -> chain ID:position:atom name
        --atom2 ATOM2         atom of the ligand to follow in this format -> chain ID:position:atom name
        --cpus CPUS           Include the number of cpus desired
        --cu                  used if there are copper in the system
        --test                Used if you want to run a test before
        --nord                used if LSF is the utility managing the jobs
        --multiple            if you want to mutate 2 residue in the same pdb
        --seed SEED           Include the seed number to make the simulation reproducible
        --dir DIR             The name of the folder for all the simulations
        --pdb_dir PDB_DIR     The name for the mutated pdb folder
        --hydrogen            leave it to default
        --consec              Consecutively mutate the PDB file for several rounds
        
