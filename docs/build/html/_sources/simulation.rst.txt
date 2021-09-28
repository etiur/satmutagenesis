Running the PELE simulations
*******************************
This script combines the functions and classes from all the other modules to create the PDB files and the yaml files, then it creates different subprocesses, as many as the number of yaml files.

.. code-block:: bash
    
    $ python -m satumut.simulation --help
    
.. code-block:: bash

    usage: simulation.py [-h] -i INPUT [-p POSITION [POSITION ...]] -lc LIGCHAIN
                     -ln LIGNAME [-at ATOMS [ATOMS ...]]
                     [-cpm CPUS_PER_MUTANT] [-tcpus TOTAL_CPUS] [-po]
                     [-fa POLARIZATION_FACTOR] [-n] [-m] [-s SEED] [-d DIR]
                     [-pd PDB_DIR] [-hy] [-co] [-st STEPS] [--dpi DPI]
                     [-tr TRAJECTORY] [--out OUT] [--plot PLOT]
                     [-an {energy,distance,both}] [--thres THRES]
                     [-sm SINGLE_MUTAGENESIS] [-PR PLURIZYME_AT_AND_RES]
                     [-r RADIUS] [-f FIXED_RESIDS [FIXED_RESIDS ...]] [-re]
                     [-are] [-x] [-cd CATALYTIC_DISTANCE]
                     [-tem TEMPLATE [TEMPLATE ...]] [-sk SKIP [SKIP ...]]
                     [-rot ROTAMERS [ROTAMERS ...]] [-e] [-l]
                     [-da DIHEDRAL_ATOMS [DIHEDRAL_ATOMS ...]] [-im {R,S}]
                     [-tu TURN] [-en ENERGY_THRESHOLD] [--QM QM]
                     [-br BOX_RADIUS]
                     [-mut {ALA,CYS,GLU,ASP,GLY,PHE,ILE,HIS,LYS,MET,LEU,ASN,GLN,PRO,SER,ARG,THR,TRP,VAL,TYR} [{ALA,CYS,GLU,ASP,GLY,PHE,ILE,HIS,LYS,MET,LEU,ASN,GLN,PRO,SER,ARG,THR,TRP,VAL,TYR} ...]]
                     [-cst {1,2}] [-pw {Binding Energy,currentEnergy}]

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
            -cpm CPUS_PER_MUTANT, --cpus_per_mutant CPUS_PER_MUTANT
                        Include the number of cpus desired
            -tcpus TOTAL_CPUS, --total_cpus TOTAL_CPUS
                        Include the number of cpus desired
            -po, --polarize_metals
                        used if there are metals in the system
            -fa POLARIZATION_FACTOR, --polarization_factor POLARIZATION_FACTOR
                        The number to divide the charges
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
            -tr TRAJECTORY, --trajectory TRAJECTORY
                        Set how many PDBs are extracted from the trajectories
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
                        and 1 letter code can be used. You can even specify
                        the protonated states
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
            -re, --restart        to place the restart flag
            -are, --adaptive_restart
                        to place the adaptive restart flag
            -x, --xtc             Change the pdb format to xtc
            -cd CATALYTIC_DISTANCE, --catalytic_distance CATALYTIC_DISTANCE
                        The distance considered to be catalytic
            -tem TEMPLATE [TEMPLATE ...], --template TEMPLATE [TEMPLATE ...]
                        Path to external forcefield templates
            -sk SKIP [SKIP ...], --skip SKIP [SKIP ...]
                        skip the processing of ligands by PlopRotTemp
            -rot ROTAMERS [ROTAMERS ...], --rotamers ROTAMERS [ROTAMERS ...]
                        Path to external rotamers templates
            -e, --equilibration   Set equilibration
            -l, --log             write logs
            -da DIHEDRAL_ATOMS [DIHEDRAL_ATOMS ...], --dihedral_atoms DIHEDRAL_ATOMS [DIHEDRAL_ATOMS ...]
                        The 4 atom necessary to calculate the dihedrals in
                        format chain id:res number:atom name
            -im {R,S}, --improve {R,S}
                        The enantiomer that should improve
            -tu TURN, --turn TURN
                        the round of plurizyme generation, not needed for the
                        1st round
            -en ENERGY_THRESHOLD, --energy_threshold ENERGY_THRESHOLD
                        The number of steps to analyse
            --QM QM               The path to the QM charges
            -br BOX_RADIUS, --box_radius BOX_RADIUS
                        Radius of the exploration box
            -mut {ALA,CYS,GLU,ASP,GLY,PHE,ILE,HIS,LYS,MET,LEU,ASN,GLN,PRO,SER,ARG,THR,TRP,VAL,TYR} , 
            --mutation {ALA,CYS,GLU,ASP,GLY,PHE,ILE,HIS,LYS,MET,LEU,ASN,GLN,PRO,SER,ARG,THR,TRP,VAL,TYR}
                        The aminoacid in 3 letter code
            -cst {1,2}, --conservative {1,2}
                        How conservative should the mutations be, choises are
                        1 and 2
            -pw {Binding Energy,currentEnergy}, --profile_with {Binding Energy,currentEnergy}
                        The metric to generate the pele profiles with

                        
| It has all the required arguments for the classes and functions from the other modules if the flags ``-sm`` and ``-PR`` are present it runs ``plurizymes_simulation`` instead of ``saturated_simulation``. 
| An example would be:

.. code-block:: bash
    
    $ python -m satumut.simulation --input PK2_F454T.pdb --position A:454 --ligchain L --ligname ANL --atoms C:1:CU L:1:N1 --cpus 5 -po --test

    
