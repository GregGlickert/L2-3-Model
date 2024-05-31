#!/bin/bash
#SBATCH --job-name=block1_short
#SBATCH --output=output_short.txt
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4

mpiexec nrniv -mpi -python batch_test2.py
