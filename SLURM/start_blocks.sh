#!/bin/bash
#SBATCH --job-name=MAIN
#SBATCH --output=start_blocks.txt
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1

python main.py