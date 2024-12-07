#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=Global_Spike
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=1000G
#SBATCH --partition=short
#SBATCH --qos=priority

module load OpenBLAS/0.3.24-GCC-13.2.0  # Load the default OpenBLAS module
module load R


# Set environment variables for multithreading
export OMP_NUM_THREADS=12         # Set OpenMP threads to 48
export OPENBLAS_NUM_THREADS=12    # Set OpenBLAS threads to 48

# Optional: Verify loaded modules and environment variables
module list
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"



Rscript Scripts/Haplotypes_Spike_BLAS.R --count.mut_threshold 8769 --Global


