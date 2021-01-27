The manual
***********

| The main usage of ``satumut`` is to mutate a given position within a protein to all the other 19 aminoacids and to facilitate their analysis and the effects on protein-ligand interactions through PELE simulations by automating the file creation and simulation launching all together. 
| As a results, it outputs 19 PDBs + 1 PDB for the wildtype and the correspoding files for the PELE simulations in marenostrum or NORD, then it launches the files in these HPCs

Introduction
===================
| After the download from the `repository <https://github.com/etiur/satumut>`_ you can readily use the pakcgae through the command line to generate the different files and lanch simulations on Marenostrum or Nord.
| Let's see the necessary arguments

.. code-block:: bash

    $ python -m satumut --help

.. code-block:: bash

    usage: main.py [-h] --input INPUT --position POSITION [POSITION ...]
               --ligchain LIGCHAIN --ligname LIGNAME --atom1 ATOM1 --atom2
               ATOM2 [--cpus CPUS] [--cu] [--test] [--nord] [--multiple]
               [--seed SEED] [--dir DIR] [--pdb_dir PDB_DIR] [--hydrogen]
               [--consec]

    Generate the mutant PDB and the corresponding running files

    optional arguments:
        -h, --help            Show this help message and exit
        --input INPUT         Include PDB file's path
        --position POSITION [POSITION ...]
                              Include one or more chain IDs and positions -> chain ID:position
        --ligchain LIGCHAIN   Include the chain ID of the ligand
        --ligname LIGNAME     The ligand residue name
        --atom1 ATOM1         Atom of the residue to follow in this format -> chain ID:position:atom name
        --atom2 ATOM2         Atom of the ligand to follow in this format -> chain ID:position:atom name
        --cpus CPUS           Include the number of cpus desired
        --cu                  Used if there are copper in the system
        --test                Used if you want to run a test before
        --nord                used if LSF is the utility managing the jobs
        --multiple            if you want to mutate 2 residue in the same pdb
        --seed SEED           Include the seed number to make the simulation reproducible
        --dir DIR             The name of the folder for all the simulations
        --pdb_dir PDB_DIR     The name for the mutated pdb folder
        --hydrogen            Leave it to default
        --consec              Consecutively mutate the PDB file for several rounds
        --sbatch              True if you want to lanch the simulation right after
                              creating the slurm file
        --steps STEPS        The number of PELE steps
        
The first 6 arguments are necessary and the rest are optional, for example:

.. code-block:: bash

    $ python -m satumut --input PK2_F454T.pdb --position A:454 --ligchain 'L' --ligname 'ANL' --atom1 "C:1:CU" --atom2 "L:1:N1" --cu --test

The code will produce a slurm file ``.sh`` and will lanch it as a job in marenostrum, then all the other files will be generated and the simulations be started by the job.
    
    
Analysis
=========
Once the simulation has been lanched, wait until the results from the simulations are generated and then you can start the analysis with the ``analysis module`` in the command line.

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
        --analyse {energy,distance,all}
                              The metric to measure the improvement of the system
        --cpus CPUS           Include the number of cpus desired
        --thres THRES         The threshold for the improvement which will affect
                              what will be included in the summary
                              
| Given a input file with the path to the folders where the PELE simulation results are stored, which is generated automatically by the main script, it will search within the       folders and generate several plots by comparing the mutations with the wildtype. 
| Then it will create a summary in **PDF format** with all the best mutations according to user defined threshold and metric of choice (energy, distance or both).

.. code-block:: bash

    $ python -m satumut.analysis --inp folder_names.txt

