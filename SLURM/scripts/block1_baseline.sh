#!/bin/bash
#SBATCH --job-name=block1_baseline
#SBATCH --output=output_baseline.txt
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4

mpiexec nrniv -mpi -python batch_test1.py
