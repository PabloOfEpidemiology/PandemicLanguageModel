#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --job-name=FG_NSP_seq_extraction
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000G
#SBATCH --partition=short
#SBATCH --qos=priority

# Define the output directory
OUTPUT_DIR="Haplotypes/24A_FG_NSP_Global_235584"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define input file paths for better readability and easier maintenance
METADATA_FILE="${OUTPUT_DIR}/NSP_Global_235584_metadata_representative_all_countries_with_date.txt"
FG_RNA_INPUT="Input_data/FG/aligned_sequences.fasta"

# Define output file paths using the OUTPUT_DIR variable
FG_RNA_OUTPUT="${OUTPUT_DIR}/NSP_Global_2355848_rna.txt"
FG_FILTERED_OUTPUT="${OUTPUT_DIR}/filtered_haplotypes.txt"

# Add RNA-level sequences of Spike to the previously generated amino acid output
python Scripts/Post_Process_Haplotypes/FG_add_rna.py "$METADATA_FILE" "$FG_RNA_INPUT" > "$FG_RNA_OUTPUT"


#Filter out all quality RNA sequences
python Scripts/Post_Process_Haplotypes/filter_low_q_RNA.py "$FG_RNA_OUTPUT" "$FG_FILTERED_OUTPUT" --FG
