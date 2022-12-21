#!/bin/bash
#SBATCH --job-name=testSimulation1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real1.out
#SBATCH --error=real1.err
time R CMD BATCH test1.R test1.out
