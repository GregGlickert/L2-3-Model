#!/bin/bash
#SBATCH --job-name=block1_baseline
#SBATCH --output=output_test/block1/baseline/%x_%j.out
#SBATCH --error=output_test/block1/baseline/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=2

export OUTPUT_DIR=output_test/block1/baseline

mpiexec nrniv -mpi -python batch_test1.py
