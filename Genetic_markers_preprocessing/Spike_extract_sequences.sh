#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --job-name=Spike_seq_extraction
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=350G
#SBATCH --partition=short
#SBATCH --qos=priority

# Define the output directory
OUTPUT_DIR="Haplotypes/24A_Spike_Countries_20272"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define input file paths for better readability and easier maintenance
METADATA_FILE="${OUTPUT_DIR}/Countries_20272_metadata_representative_all_countries_with_date.txt"
SPIKE_AA_INPUT="Input_data/Spike_only/gene_S.translation.fasta"
SPIKE_RNA_INPUT="Input_data/raw_GISAID_data/spikenuc1026.fasta"

# Define output file paths using the OUTPUT_DIR variable
SPIKE_AA_OUTPUT="${OUTPUT_DIR}/Spike_Countries_20272_Spike_aa.txt"
SPIKE_AA_RNA_OUTPUT="${OUTPUT_DIR}/Spike_Countries_20272_Spike_aa_rna.txt"
SPIKE_FILTERED_OUTPUT="${OUTPUT_DIR}/filtered_haplotypes.txt"

# Add amino acid-level sequences of Spike to haplotype metadata
python Scripts/Post_Process_Haplotypes/Spike_add_aa.py "$METADATA_FILE" "$SPIKE_AA_INPUT" > "$SPIKE_AA_OUTPUT"

# Add RNA-level sequences of Spike to the previously generated amino acid output
python Scripts/Post_Process_Haplotypes/Spike_add_rna.py "$SPIKE_AA_OUTPUT" "$SPIKE_RNA_INPUT" > "$SPIKE_AA_RNA_OUTPUT"


#Filter out all quality RNA sequences
python Scripts/Post_Process_Haplotypes/filter_low_q_RNA.py "$SPIKE_AA_RNA_OUTPUT" "$SPIKE_FILTERED_OUTPUT" --Spike
