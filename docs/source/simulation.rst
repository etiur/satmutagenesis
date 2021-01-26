Running the PELE simulations
*******************************
This script combines the functions and classes from all the other modules to create the PDB files and the yaml files, then it creates different subprocesses, as many as the number of yaml files.

.. code-block:: bash
    
    $ python -m satumut.simulation --help
    
.. code-block:: bash

    usage: simulation.py [-h] --input INPUT --position POSITION [POSITION ...]
                     --ligchain LIGCHAIN --ligname LIGNAME --atom1 ATOM1
                     --atom2 ATOM2 [--cpus CPUS] [--cu] [--test] [--nord]
                     [--multiple] [--seed SEED] [--dir DIR]
                     [--pdb_dir PDB_DIR] [--hydrogen] [--consec]

    Generate the mutant PDB and the corresponding running files

    optional arguments:
        -h, --help            show this help message and exit
        --input INPUT         Include PDB file's path
        --position POSITION [POSITION ...]
                              Include one or more chain IDs and positions -> Chain
                              ID:position
        --ligchain LIGCHAIN   Include the chain ID of the ligand
        --ligname LIGNAME     The ligand residue name
        --atom1 ATOM1         atom of the residue to follow in this format -> chain
                              ID:position:atom name
        --atom2 ATOM2         atom of the ligand to follow in this format -> chain
                              ID:position:atom name
        --cpus CPUS           Include the number of cpus desired
        --cu                  used if there are copper in the system
        --test                Used if you want to run a test before
        --nord                used if LSF is the utility managing the jobs
        --multiple            if you want to mutate 2 residue in the same pdb
        --seed SEED           Include the seed number to make the simulation
                              reproducible
        --dir DIR             The name of the folder for all the simulations
        --pdb_dir PDB_DIR     The name for the mutated pdb folder
        --hydrogen            leave it to default
        --consec              Consecutively mutate the PDB file for several rounds

| It has all the required arguments for the classes and functions from the other modules. 
| An example would be:

.. code-block:: bash
    
    $ python -m satumut.simulation --input PK2_F454T.pdb --position A:454 --ligchain L --ligname ANL --atom1 C:1:CU --atom2 L:1:N1 --cpus 5 --cu --test

    
