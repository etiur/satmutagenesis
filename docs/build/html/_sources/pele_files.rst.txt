The creation of files for PELE simulations
********************************************
The second script is used when you have the PDB files and you want to generate the control files for lanching the PELE simulation for each of the PDBs.

.. code-block:: bash

    $ python -m satumut.mutate_pdb --help

.. code-block:: bash

    usage: pele_files.py [-h] --folder FOLDER -lc LIGCHAIN -ln LIGNAME -at ATOMS
                 [ATOMS ...] [--cpus CPUS] [-po] [-fa POLARIZATION_FACTOR]
                 [-t] [-n] [-s SEED] [-st STEPS]

    Generate running files for PELE

    optional arguments:
        -h, --help            show this help message and exit
        --folder FOLDER       An iterable of the path to different pdb files, a name
                    of the folder or a file with the path to the different
                    pdb files
        -lc LIGCHAIN, --ligchain LIGCHAIN
                    Include the chain ID of the ligand
        -ln LIGNAME, --ligname LIGNAME
                    The ligand residue name
        -at1 ATOM1, --atom1 ATOM1
                        atom of the residue to follow in this format -> chain
                        ID:position:atom name
        -at2 ATOM2, --atom2 ATOM2
                        atom of the ligand to follow in this format -> chain
                        ID:position:atom name
        --cpus CPUS           Include the number of cpus desired
        -po, --polarize_metals
                        used if there are metals in the system
        -fa POLARIZATION_FACTOR, --polarization_factor POLARIZATION_FACTOR
                          The number to divide the charges
        -t, --test            Used if you want to run a test before
        -n, --nord            used if LSF is the utility managing the jobs
        -s SEED, --seed SEED  Include the seed number to make the simulation
                    reproducible
        -st STEPS, --steps STEPS
                    The number of PELE steps
                             
The ``--folder`` flag accepts folders where the PDB files are located, files containing the path to the different PDBs (one path per line) or lists of paths to the PDBs, for example:

.. code-block:: bash

    $ python -m satumut.mutate_pdb --folder pdb_files --ligchain 'L' --ligname 'ANL' --atom1 "C:1:CU" --atom2 "L:1:N1" -po --test
    
As a result, it will create the 2 files necessary to launch a PELE simulation in marenostrum or nord for each of the PDBs in the folder.

