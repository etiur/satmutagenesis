Running the PELE simulations
*******************************
This script combines the functions and classes from all the other modules to create the PDB files and the yaml files, then it creates different subprocesses, as many as the number of yaml files.

.. code-block:: bash
    
    $ python -m satumut.simulation --help
    
.. code-block:: bash

    usage: simulation.py [-h] -i INPUT [-p POSITION [POSITION ...]] -lc LIGCHAIN
                     -ln LIGNAME -at ATOMS [ATOMS ...] [--cpus CPUS] [-po]
                     [-fa POLARIZATION_FACTOR] [-t] [-n] [-m] [-s SEED] [-d DIR] 
                     [-pd PDB_DIR] [-hy]
                     [-co]   [-st STEPS] [--dpi DPI] [--box BOX]
                     [--traj TRAJ] [--out OUT] [--plot PLOT]
                     [-an {energy,distance,both}] [--thres THRES]
                     [-sm SINGLE_MUTAGENESIS] [-PR PLURIZYME_AT_AND_RES]
                     [-r RADIUS] [-f FIXED_RESIDS [FIXED_RESIDS ...]] [-pa]

    Generate the mutant PDB and the corresponding running files

    optional arguments:
        -h, --help            show this help message and exit
        -i INPUT, --input INPUT
                        Include PDB file's path
        -p POSITION [POSITION ...], --position POSITION [POSITION ...]
                        Include one or more chain IDs and positions -> Chain
                        ID:position
        -lc LIGCHAIN, --ligchain LIGCHAIN
                        Include the chain ID of the ligand
        -ln LIGNAME, --ligname LIGNAME
                        The ligand residue name
        -at ATOMS [ATOMS ...], --atoms ATOMS [ATOMS ...]
                        Series of atoms of the residues to follow in this
                        format -> chain ID:position:atom name
        --cpus CPUS           Include the number of cpus desired
        -po, --polarize_metals
                            used if there are metals in the system
        -fa POLARIZATION_FACTOR, --polarization_factor POLARIZATION_FACTOR
                              The number to divide the charges
        -t, --test            Used if you want to run a test before
        -n, --nord            used if LSF is the utility managing the jobs
        -m, --multiple        if you want to mutate 2 residue in the same pdb
        -s SEED, --seed SEED  Include the seed number to make the simulation
                        reproducible
        -d DIR, --dir DIR     The name of the folder for all the simulations
        -pd PDB_DIR, --pdb_dir PDB_DIR
                        The name for the mutated pdb folder
        -hy, --hydrogen       leave it to default
        -co, --consec         Consecutively mutate the PDB file for several rounds

        -st STEPS, --steps STEPS
                        The number of PELE steps
        --dpi DPI             Set the quality of the plots
        --box BOX             Set how many data points are used for the boxplot
        --traj TRAJ           Set how many PDBs are extracted from the trajectories
        --out OUT             Name of the summary file created at the end of the
                        analysis
        --plot PLOT           Path of the plots folder
        -an {energy,distance,both}, --analyse {energy,distance,both}
                        The metric to measure the improvement of the system
        --thres THRES         The threshold for the improvement which will affect
                        what will be included in the summary
        -sm SINGLE_MUTAGENESIS, --single_mutagenesis SINGLE_MUTAGENESIS
                        Specifiy the name of the residue that you want the
                        original residue to be mutated to. Both 3 letter code
                        and 1 letter code can be used. You can even specify the protonated states
        -PR PLURIZYME_AT_AND_RES, --plurizyme_at_and_res PLURIZYME_AT_AND_RES
                        Specify the chain ID, residue number and the PDB atom
                        name thatwill set the list of the neighbouring
                        residues for thenext round. Example: chain
                        ID:position:atom name
        -r RADIUS, --radius RADIUS
                        The radius around the selected atom to search for the
                        other residues
        -f FIXED_RESIDS [FIXED_RESIDS ...], --fixed_resids FIXED_RESIDS [FIXED_RESIDS ...]
                        Specify the list of residues that you don't wantto
                        have mutated (Must write the list of residuenumbers)
        -pa, --pele_analysis  if you want to turn on the analysis by PELE
                        
| It has all the required arguments for the classes and functions from the other modules if the falgs ``-sm`` and ``-PR`` are present it runs ``plurizymes_simulation`` instead of ``saturated_simulation``. 
| An example would be:

.. code-block:: bash
    
    $ python -m satumut.simulation --input PK2_F454T.pdb --position A:454 --ligchain L --ligname ANL --atoms C:1:CU L:1:N1 --cpus 5 -po --test

    
