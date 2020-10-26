#!/bin/bash
#SBATCH -J PELE
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --qos=debug
#SBATCH --ntasks=5

module purge
export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"
export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"
export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH
module load impi
module load intel mkl impi gcc # 2> /dev/null
module load boost/1.64.0
/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8 -m pele_platform.main ../yaml_files/GLN_38.yaml
