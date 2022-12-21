#!/bin/bash
#SBATCH --job-name=base0_1.2_400
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=real.out
#SBATCH --error=real.err
time R CMD BATCH phynotypepppower_betaf0=0_betagc=1.2_sample400.R apppower_betaf0=0_betagc=1.2_sample400.out
