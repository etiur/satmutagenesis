The creation of yaml files for PELE platform
**********************************************
The second script is used when you have the PDB files and you want to generate the yamls files for each of the PDBs.

.. code-block:: bash

    $ python -m satumut.pele_files --help

.. code-block:: bash

    usage: pele_files.py [-h] --folder FOLDER -lc LIGCHAIN -ln LIGNAME
                     [-at ATOMS [ATOMS ...]] [-cpm CPUS_PER_MUTANT]
                     [-tcpus TOTAL_CPUS] [-po] [-fa POLARIZATION_FACTOR] [-n]
                     [-s SEED] [-st STEPS] [-x] [-e]
                     [-tem TEMPLATE [TEMPLATE ...]]
                     [-rot ROTAMERS [ROTAMERS ...]] [-sk SKIP [SKIP ...]] [-l]
                     [-co] [-tu TURN] [--QM QM] [-br BOX_RADIUS]

    Generate running files for PELE

    optional arguments:
        -h, --help            show this help message and exit
        --folder FOLDER       An iterable of the path to different pdb files, a name
                        of the folder with the pdbs
        -lc LIGCHAIN, --ligchain LIGCHAIN
                        Include the chain ID of the ligand
        -ln LIGNAME, --ligname LIGNAME
                        The ligand residue name
        -at ATOMS [ATOMS ...], --atoms ATOMS [ATOMS ...]
                        Series of atoms of the residues to follow in this
                        format -> chain ID:position:atom name
        -cpm CPUS_PER_MUTANT, --cpus_per_mutant CPUS_PER_MUTANT
                        Include the number of cpus desired
        -tcpus TOTAL_CPUS, --total_cpus TOTAL_CPUS
                        Include the number of cpus desired
        -po, --polarize_metals
                        used if there are metals in the system
        -fa POLARIZATION_FACTOR, --polarization_factor POLARIZATION_FACTOR
                        The number to divide the charges
        -n, --nord            used if LSF is the utility managing the jobs
        -s SEED, --seed SEED  Include the seed number to make the simulation
                        reproducible
        -st STEPS, --steps STEPS
                        The number of PELE steps
        -x, --xtc             Change the pdb format to xtc
        -e, --equilibration   Set equilibration
        -tem TEMPLATE [TEMPLATE ...], --template TEMPLATE [TEMPLATE ...]
                        Path to external forcefield templates
        -rot ROTAMERS [ROTAMERS ...], --rotamers ROTAMERS [ROTAMERS ...]
                        Path to external rotamers templates
        -sk SKIP [SKIP ...], --skip SKIP [SKIP ...]
                        skip the processing of ligands by PlopRotTemp
        -l, --log             write logs
        -co, --consec         Consecutively mutate the PDB file for several rounds
        -tu TURN, --turn TURN
                        the round of plurizyme generation, not needed for the
                        1st round
        --QM QM               The path to the QM charges
        -br BOX_RADIUS, --box_radius BOX_RADIUS
                        Radius of the exploration box

The ``-br`` flag modifies the radius of the box where induce_fir simulation happens

.. code-block:: bash

    $ python -m satumut.pele_files --folder pdb_files --ligchain 'L' --ligname 'ANL' --atoms "C:1:CU" "L:1:N1" -po --test
    
As a result, it will create the yaml files which are the input files for the PELE platform. 

