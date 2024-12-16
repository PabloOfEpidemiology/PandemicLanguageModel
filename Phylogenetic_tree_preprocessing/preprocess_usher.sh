#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=usher
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000G
#SBATCH --partition=short
#SBATCH --qos=priority

#python generate_epiToPublicAndDate.py
python scripts/multiCPU_preprocess_usher.py \
    --tree-file-in results/gisaidAndPublic.2024-12-10.masked.pb \
    --gisaid-metadata-file-in results/metadata.tsv.gz \
    --usher-metadata-file-in results/gisaidAndPublic.2024-12-10.metadata.tsv \
    --max-num-clades 12000,20000  \
    --pango-loss 0,100000 \
    --skip-features \
    --num-cpus 8

