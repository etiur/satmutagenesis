The generation of mutations
****************************
Now let's look at the individual scripts starting with the mutate_pdb.py which can be used if you only want to generate the PDB files of the mutations. It outputs 19 PDBs + 1 PDB for the wildtype.

.. code-block:: bash

    $ python -m saturated_mutagenesis.mutate_pdb --help
    
.. code-block:: bash

    usage: mutate_pdb.py [-h] -i INPUT -p POSITION [POSITION ...] [-m] [-hy] [-co]
                     [-pd PDB_DIR] [-sm SINGLE_MUTAGENESIS]

    Performs saturated mutagenesis given a PDB file

    optional arguments:
        -h, --help            show this help message and exit
        -i INPUT, --input INPUT
                            Include PDB file's path
        -p POSITION [POSITION ...], --position POSITION [POSITION ...]
                            Include one or more chain IDs and positions -> Chain
                            ID:position
        -m, --multiple        if you want to mutate 2 residue in the same pdb
        -hy, --hydrogen       leave it to default
        -co, --consec         Consecutively mutate the PDB file for several rounds
        -pd PDB_DIR, --pdb_dir PDB_DIR
                            The name for the mutated pdb folder
        -sm SINGLE_MUTAGENESIS, --single_mutagenesis SINGLE_MUTAGENESIS
                            Specifiy the name of the residue that you want the
                            original residue to be mutated to. Both 3 letter code
                            and 1 letter code can be used.
        
There are 2 necessary arguments, input or the wildtype PDB file and the position or positions to mutate, for example:

.. code-block:: bash

    $  python -m saturated_mutagenesis.mutate_pdb --input PK2_F454T.pdb --position A:454 A:135

| If the flag ``--multiple`` is present and there are 2 positions to mutate to, it will generate 400 PDBs as a result of the combination of 20 X 20 mutations.
| If the flag ``--consec`` is present, it means that you have mutated 20 residues the first time, and you are selecting 1 of those 20 mutations to start the second round fo mutations. The flag preserves the name of the PDB file so you know from which PDB file it came from the second mutation.

