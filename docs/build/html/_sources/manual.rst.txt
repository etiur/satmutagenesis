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

    usage: __main__.py [-h] -i INPUT -p POSITION [POSITION ...] -lc LIGCHAIN -ln
                   LIGNAME -at1 ATOM1 -at2 ATOM2 [--cpus CPUS] [-po]
                   [-fa POLARIZATION_FACTOR] [-t] [-n] [-m] [-s SEED] [-d DIR]
                   [-pd PDB_DIR] [-hy] [-co] [-st STEPS]

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
        
The first 6 arguments are necessary and the rest are optional, for example:

.. code-block:: bash

    $ python -m satumut --input PK2_F454T.pdb --position A:454 --ligchain 'L' --ligname 'ANL' --atom1 "C:1:CU" --atom2 "L:1:N1" -po --test


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

