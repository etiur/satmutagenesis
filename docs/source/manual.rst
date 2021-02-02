The manual
***********

| The main usage of ``satumut`` is to mutate a given position within a protein to all the other 19 aminoacids and to facilitate their analysis and the effects on protein-ligand interactions through PELE simulations by automating the file creation and simulation launching all together. 
| As a results, it outputs 19 PDBs + 1 PDB for the wildtype and the correspoding files for the PELE simulations in marenostrum or NORD, then it launches the files in these HPCs

Introduction
===================
| After the download from the `repository <https://github.com/etiur/satumut>`_ you can readily use the pakcgae through the command line to generate the different files, lanch the simulations on Marenostrum or Nord and analyse the simulations through boxplots, scatter plots and trajectories.
| Let's see the necessary arguments

.. code-block:: bash

    $ python -m satumut --help

.. code-block:: bash

    usage: __main__.py [-h] -i INPUT [-p POSITION [POSITION ...]] -lc LIGCHAIN -ln
                   LIGNAME -at ATOMS [ATOMS ...] [--cpus CPUS] [--cu] [-t]
                   [-n] [-m] [-s SEED] [-d DIR] [-pd PDB_DIR] [-hy] [-co]
                   [-sb] [-st STEPS] [--dpi DPI] [--box BOX] [--traj TRAJ]
                   [--out OUT] [--plot PLOT] [-an {energy,distance,both}]
                   [--thres THRES] [-sm SINGLE_MUTAGENESIS]
                   [-PR PLURIZYME_AT_AND_RES] [-r RADIUS]
                   [-f FIXED_RESIDS [FIXED_RESIDS ...]]

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
        --cu                  used if there are copper in the system
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
        -sb, --sbatch         True if you want to lanch the simulation right after
                              creating the slurm file
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
                              and 1 letter code can be used.
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
                              have mutated (Must write the list of residue
                              positionnumbers)
                              
The first 6 arguments are necessary and the rest are optional, for example:

.. code-block:: bash

    $ python -m satumut --input PK2_F454T.pdb --position A:454 --ligchain 'L' --ligname 'ANL' --atoms "C:1:CU" "L:1:N1" --cu --test

The code will produce a slurm file ``.sh`` and will lanch it as a job in marenostrum, then all the other files will be generated and the simulations be started by the job.
    
Analysis
=========
Once the simulation has been lanched, The job will wait until the results from the simulations are generated and then it will start with the analysis with the ``analysis module``.

.. code-block:: bash

    $ python -m satumut.analysis --help
    
.. code-block:: bash

    usage: analysis.py [-h] --inp INP [--dpi DPI] [--box BOX] [--traj TRAJ]
                   [--out OUT] [--folder FOLDER]
                   [--analyse {energy,distance,all}] [--cpus CPUS]
                   [--thres THRES]

    Analyse the different PELE simulations and create plots

    optional arguments:
        -h, --help            Show this help message and exit
        --inp INP             Include a file or list with the path to the folders
                              with PELE simulations inside
        --dpi DPI             Set the quality of the plots
        --box BOX             Set how many data points are used for the boxplot
        --traj TRAJ           Set how many PDBs are extracted from the trajectories
        --out OUT             Name of the summary file created at the end of the
                              analysis
        --folder FOLDER       Name of the plots folder
        --analyse {energy,distance,both}
                              The metric to measure the improvement of the system
        --cpus CPUS           Include the number of cpus desired
        --thres THRES         The threshold for the improvement which will affect
                              what will be included in the summary
                              
| Given a input file with the path to the folders where the PELE simulation results are stored, which is generated automatically by the main script, it will search within the       folders and generate several plots by comparing the mutations with the wildtype. 
| Then it will create a summary in **PDF format** with all the best mutations according to user defined threshold and metric of choice (energy, distance or both).

.. code-block:: bash

    $ python -m satumut.analysis --inp folder_names.txt

