#!/bin/bash

#SBATCH --partition=amd
##SBATCH -w SG12
##SBATCH --exclude=SG14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
##SBATCH --job-name=raspa3
#SBATCH --export=ALL

cd $SLURM_SUBMIT_DIR

module purge
module add raspa3-mbx/E_print

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

raspa3-mbx


