import sys
import csv
import os
from Bio import SeqIO as sio

def to_dict(path=None, file_format='fasta', unique=False):
    if path is None:
        sys.stderr.write('Path to FASTA file required\n')
        return {}
    res = {}
    if unique:
        sys.stderr.write('Removed duplicate name occurrences without checking sequences\n')
        names = set()
        for i, record in enumerate(sio.parse(path, file_format), 1):
            if record.description not in names:
                names.add(record.description)
                res[record.description] = str(record.seq)
            if i % 100000 == 0:
                sys.stderr.write(f'Parsed {i} sequences from FASTA\n')
    else:
        for i, record in enumerate(sio.parse(path, file_format), 1):
            res[record.description] = str(record.seq)
            if i % 100000 == 0:
                sys.stderr.write(f'Parsed {i} sequences from FASTA\n')
    return res

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python Scripts/add_Spike_aa.py metadata.(csv|txt) sequences.fasta > output.(csv|txt)\n")
        sys.exit(1)
    
    metadata_path = sys.argv[1]
    fasta_path = sys.argv[2]
    
    # Determine delimiter based on metadata file extension
    _, ext = os.path.splitext(metadata_path)
    ext = ext.lower()
    if ext == '.csv':
        delimiter = ','
    elif ext == '.txt':
        delimiter = '\t'
    else:
        sys.stderr.write("Metadata file must be .csv or .txt\n")
        sys.exit(1)
    
    # Create a dictionary of sequences
    seq_d = to_dict(path=fasta_path, unique=True)
    
    not_found_count = 0  # Initialize counter for sequences not found
    
    with open(metadata_path, newline='') as metadata_f:
        reader = csv.DictReader(metadata_f, delimiter=delimiter)
        fieldnames = reader.fieldnames + ['aa']  # Add new column for Spike amino acid sequences
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        
        for i, row in enumerate(reader, 1):
            name = row.get('seqName', '')
            sequence = seq_d.get(name)
            if sequence:
                row['aa'] = sequence
            else:
                row['aa'] = "Sequence_Not_Found"
                not_found_count += 1
            writer.writerow(row)
            if i % 1000 == 0:
                sys.stderr.write(f'Processed {i} rows\n')
    
    sys.stderr.write(f'Total processed rows: {i}\n')
    sys.stderr.write(f'Total sequences not found: {not_found_count}\n')

if __name__ == "__main__":
    main()

