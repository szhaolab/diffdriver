#!/bin/bash
#SBATCH --job-name=dd0_0_1200
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH phynotypepppower_betaf0=0_betagc=0_sample1200.R apppower_betaf0=0_betagc=0_sample1200.out
