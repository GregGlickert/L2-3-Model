#!/bin/bash
#SBATCH --job-name=block2_long
#SBATCH --output=output_test/block2/long/%x_%j.out
#SBATCH --error=output_test/block2/long/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=2

export OUTPUT_DIR=output_test/block2/long

mpiexec nrniv -mpi -python batch_test3.py