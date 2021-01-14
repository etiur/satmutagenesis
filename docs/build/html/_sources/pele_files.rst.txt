The creation of files for PELE simulations
********************************************
The second script is used when you have the PDB files and you want to generate the control files for lanching the PELE simulation for each of the PDBs.

.. code-block:: bash

    $ python -m saturated_mutagenesis.mutate_pdb --help

.. code-block:: bash

    usage: pele_files.py [-h] --folder FOLDER --ligchain LIGCHAIN --ligname
                     LIGNAME --atom1 ATOM1 --atom2 ATOM2 [--cpus CPUS] [--cu]
                     [--test] [--nord] [--seed SEED]

    Generate running files for PELE

    optional arguments:
        -h, --help           show this help message and exit
        --folder FOLDER      An iterable of the path to different pdb files, a name
                             of the folder or a file of the path to the different
                             pdb files
        --ligchain LIGCHAIN  Include the chain ID of the ligand
        --ligname LIGNAME    The ligand residue name
        --atom1 ATOM1        atom of the residue to follow in this format -> chain
                             ID:position:atom name
        --atom2 ATOM2        atom of the ligand to follow in this format -> chain
                             ID:position:atom name
        --cpus CPUS          Include the number of cpus desired
        --cu                 used if there are copper in the system
        --test               Used if you want to run a test before
        --nord               used if LSF is the utility managing the jobs
        --seed SEED          Include the seed number to make the simulation
                             reproducible
                             
The ``--folder`` flag accepts folders where the PDB files are located, files containing the path to the different PDBs (one path per line) or lists of paths to the PDBs, for example:

.. code-block:: bash

    $ python -m saturated_mutagenesis.mutate_pdb --folder pdb_files --ligchain 'L' --ligname 'ANL' --atom1 "C:1:CU" --atom2 "L:1:N1" --cu --test
    
As a result, it will create the 2 files necessary to launch a PELE simulation in marenostrum or nord for each of the PDBs in the folder.

