#!/bin/bash
#SBATCH --job-name=testSimulation2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real2.out
#SBATCH --error=real2.err
time R CMD BATCH test2.R test2.out
