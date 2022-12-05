#!/bin/bash
#SBATCH --job-name=dd1_0_600
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH phynotypepppower_betaf0=1_betagc=0_sample600.R apppower_betaf0=1_betagc=0_sample600.out
