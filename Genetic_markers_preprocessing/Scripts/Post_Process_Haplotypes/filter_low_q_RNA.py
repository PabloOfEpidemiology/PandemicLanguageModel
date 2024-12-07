#!/usr/bin/env python3

import sys
import csv
import argparse

def calculate_length(sequence):
    """
    Calculates the length of the RNA sequence.
    """
    return len(sequence)

def calculate_low_quality_fraction(sequence):
    """
    Calculates the fraction of low-quality nucleotides in the RNA sequence.
    Low-quality nucleotides are those not in the set {'A', 'T', 'G', 'C', 'U', '-'}.
    """
    valid_nucleotides = set(['A', 'T', 'G', 'C', 'U', '-'])
    total_length = len(sequence)
    if total_length == 0:
        return 0.0

    low_quality_count = sum(1 for nucleotide in sequence.upper() if nucleotide not in valid_nucleotides)
    fraction = low_quality_count / total_length
    return fraction  # Return as float

def main():
    parser = argparse.ArgumentParser(description='Filter RNA sequences based on quality control criteria.')
    parser.add_argument('input_file', help='Input file (.txt or .csv)')
    parser.add_argument('output_file', help='Output file (.txt or .csv)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--FG', action='store_true', help='Full genome sequences')
    group.add_argument('--Spike', action='store_true', help='Spike sequences')
    group.add_argument('--Multigene', action='store_true', help='Multigene sequences')
    args = parser.parse_args()

    # Set quality control thresholds based on the flag
    if args.FG:
        min_length = 28000
    elif args.Spike:
        min_length = 3500
    elif args.Multigene:
        min_length = 6635

    # Initialize counters
    total_rows = 0
    failed_qc = 0

    # Open input and output files
    with open(args.input_file, 'r', newline='') as infile, open(args.output_file, 'w', newline='') as outfile:
        # Detect delimiter (assume tab-delimited)
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for i, row in enumerate(reader, 1):
            total_rows +=1

            if args.FG or args.Spike:
                # Get the RNA sequence from 'RNA' column
                rna_sequence = row.get('RNA', '')
            elif args.Multigene:
                # Concatenate all RNA sequences from columns starting with 'RNA'
                rna_sequence = ''
                for key in row:
                    if key.startswith('RNA'):
                        rna_sequence += row[key]

            # Calculate length and low-quality fraction
            seq_length = calculate_length(rna_sequence)
            low_quality_fraction = calculate_low_quality_fraction(rna_sequence)

            # Check quality criteria
            if seq_length < min_length or low_quality_fraction > 0.05:
                failed_qc +=1
                continue  # Skip writing this row to output file

            # Write the row to output file
            writer.writerow(row)

            if i % 1000 == 0:
                sys.stderr.write(f'Processed {i} rows\n')

        sys.stderr.write(f'Total processed rows: {total_rows}\n')
        sys.stderr.write(f'Rows failed QC: {failed_qc}\n')

if __name__ == "__main__":
    main()

