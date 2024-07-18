#!/bin/sh
#SBATCH -J M1_sim 
#SBATCH -o  M1_sim.out
#SBATCH -e  M1_sim.error
#SBATCH -t 0-48:00:00  # days-hours:minutes

#SBATCH -N 1
#SBATCH -n 50 # used for MPI codes, otherwise leave at '1'
##SBATCH --ntasks-per-node=1  # don't trust SLURM to divide the cores evenly
##SBATCH --cpus-per-task=1  # cores per task; set to one if using MPI
##SBATCH --exclusive  # using MPI with 90+% of the cores you should go exclusive
#SBATCH --mem-per-cpu=4G  # memory per core; default is 1GB/core

START=$(date)
echo "Started running at $START."

export HDF5_USE_FILE_LOCKING=FALSE
unset DISPLAY

export OUTPUT_DIR=../Run-Storage/syn-weight-check/block1
mpirun nrniv -mpi -python run_network.py simulation_config_baseline.json False # args: config file, whether use coreneuron

END=$(date)
echo "Done running simulation at $END"

