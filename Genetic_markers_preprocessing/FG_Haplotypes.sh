#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=RNA_Hap
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=1000G
#SBATCH --partition=short


module load OpenBLAS/0.3.24-GCC-13.2.0  # Load the default OpenBLAS module
module load R


# Set environment variables for multithreading
export OMP_NUM_THREADS=48         # Set OpenMP threads to 48
export OPENBLAS_NUM_THREADS=48    # Set OpenBLAS threads to 48

# Optional: Verify loaded modules and environment variables
module list
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"


Rscript Scripts/Haplotypes_FG.R --count.mut_threshold 2012738 --Global --RNA


