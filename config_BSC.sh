#!/bin/bash
module purge 
module load gcc/4.9.4 ANACONDA/2019.10
export CC=gcc
conda activate /gpfs/projects/bsc72/conda_envs/saturated
