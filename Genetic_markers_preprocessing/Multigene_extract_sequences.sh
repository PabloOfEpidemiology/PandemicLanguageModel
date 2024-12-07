#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --job-name=Multigene_seq_extraction
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000G
#SBATCH --partition=short
#SBATCH --qos=priority

# Load necessary modules (if BioPython is installed via modules)
# module load python/3.8

# Define the output directory
OUTPUT_DIR="Haplotypes/24A_Multigene_NSP_Global_80072"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define input file paths for better readability and easier maintenance
METADATA_FILE="${OUTPUT_DIR}/NSP_Global_80072_metadata_representative_all_countries_with_date.txt"
aa_FASTA_INPUT_DIR="Input_data/FG/Multigene_S_ORF3a_ORF7a_ORF8_N"
RNA_FASTA_INPUT="Input_data/FG/Multigene_S_ORF3a_ORF7a_ORF8_N/aligned_sequences.fasta"


# Define output file paths using the OUTPUT_DIR variable
Multigene_aa_OUTPUT="${OUTPUT_DIR}/NSP_Global_80072_aa.txt"
Multigene_aa_RNA_OUTPUT="${OUTPUT_DIR}/NSP_Global_80072_aa_RNA.txt"
Multigene_FILTERED_OUTPUT="${OUTPUT_DIR}/filtered_haplotypes.txt"

# Add aa-level sequences of multiple genes to the metadata
python Scripts/Post_Process_Haplotypes/Multigene_add_aa.py "$METADATA_FILE" "$aa_FASTA_INPUT_DIR" > "$Multigene_aa_OUTPUT"

# Check if the previous command was successful
if [ $? -ne 0 ]; then
    echo "Multigene_add_aa.py failed. Exiting."
    exit 1
fi



# Add RNA-level sequences of multiple genes to the metadata (same logic of sequence adding as in FG, so using the same script)
python Scripts/Post_Process_Haplotypes/FG_add_rna.py "$Multigene_aa_OUTPUT" "$RNA_FASTA_INPUT" > "$Multigene_aa_RNA_OUTPUT"

# Check if the previous command was successful
if [ $? -ne 0 ]; then
    echo "FG_add_RNA.py failed. Exiting."
    exit 1
fi




# Filter out all low-quality RNA sequences
python Scripts/Post_Process_Haplotypes/filter_low_q_RNA.py "$Multigene_aa_RNA_OUTPUT" "$Multigene_FILTERED_OUTPUT" --Multigene

# Check if the filtering was successful
if [ $? -ne 0 ]; then
    echo "filter_low_q_RNA.py failed. Exiting."
    exit 1
fi

echo "Sequence extraction and filtering completed successfully."

