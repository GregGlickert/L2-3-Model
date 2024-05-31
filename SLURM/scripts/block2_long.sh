#!/bin/bash
#SBATCH --job-name=block2_long
#SBATCH --output=output_long.txt
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4

mpiexec nrniv -mpi -python batch_test3.py
