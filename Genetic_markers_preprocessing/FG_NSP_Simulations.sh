#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=parNSP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3000G
#SBATCH --partition=short
#SBATCH --qos=priority


module load R

Rscript Scripts/Simulations_FG_parallel.R --num.cores 9 --count.mut 500,1000,1500,3000,7000,10000,20000,40000,70000,100000,140000,200000,400000,600000,800000,1000000,1250000,1500000 --global --NSP


