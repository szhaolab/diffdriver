#!/bin/bash
#SBATCH --job-name=other1_0_1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH phynotypepppower_betaf0=1_betagc=0_sample1000.R apppower_betaf0=1_betagc=0_sample1000.out
