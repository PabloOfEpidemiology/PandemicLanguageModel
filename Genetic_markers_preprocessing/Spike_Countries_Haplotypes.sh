#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=22B_C
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=350G
#SBATCH --partition=short
#SBATCH --qos=priority

module load R

Rscript Scripts/Haplotypes_Spike_Countries.R --count.mut_threshold 20272 --Countries

