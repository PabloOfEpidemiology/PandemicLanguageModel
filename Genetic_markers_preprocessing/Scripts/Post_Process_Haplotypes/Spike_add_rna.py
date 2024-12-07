import sys
import csv
import os
from Bio import SeqIO as sio

def extract_accession_id(description):
    """
    Extracts the Accession.ID (e.g., EPI_ISL_402124) from the FASTA header description.
    
    Parameters:
        description (str): The description line from the FASTA header.
    
    Returns:
        str or None: The extracted Accession.ID or None if not found.
    """
    parts = description.split('|')
    for part in parts:
        if part.startswith('EPI_ISL_'):
            return part
    return None  # Return None if Accession.ID is not found

def to_dict(path=None, file_format='fasta', unique=False):
    """
    Converts a FASTA file into a dictionary with Accession.ID as keys and sequences as values.
    
    Parameters:
        path (str): Path to the FASTA file.
        file_format (str): Format of the FASTA file (default is 'fasta').
        unique (bool): If True, only the first occurrence of each Accession.ID is stored.
    
    Returns:
        dict: A dictionary mapping Accession.ID to sequences.
    """
    if path is None:
        sys.stderr.write('Path to FASTA file required\n')
        return {}
    
    res = {}
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        if unique:
            sys.stderr.write('Removed duplicate Accession.ID occurrences by keeping the first one\n')
            accession_ids = set()
            for i, record in enumerate(sio.parse(handle, file_format), 1):
                accession_id = extract_accession_id(record.description)
                if accession_id and accession_id not in accession_ids:
                    accession_ids.add(accession_id)
                    res[accession_id] = str(record.seq)
                if i % 100000 == 0:
                    sys.stderr.write(f'Parsed {i} sequences from FASTA\n')
        else:
            for i, record in enumerate(sio.parse(handle, file_format), 1):
                accession_id = extract_accession_id(record.description)
                if accession_id:
                    res[accession_id] = str(record.seq)
                if i % 100000 == 0:
                    sys.stderr.write(f'Parsed {i} sequences from FASTA\n')
    return res

def main():
    """
    Main function to process the metadata file and add RNA sequences based on Accession.ID.
    """
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python Scripts/add_rna.py metadata.(csv|txt) sequences.fasta > output.(csv|txt)\n")
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
    
    # Create a dictionary of sequences using Accession.ID as the key
    seq_d = to_dict(path=fasta_path, unique=True)
    
    not_found_count = 0  # Initialize counter for sequences not found
    
    with open(metadata_path, newline='') as metadata_f:
        reader = csv.DictReader(metadata_f, delimiter=delimiter)
        # Add the new column 'Spike_RNA' to the existing fieldnames
        fieldnames = reader.fieldnames + ['RNA']
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        
        for i, row in enumerate(reader, 1):
            accession_id = row.get('Accession.ID', '')  # Use 'Accession.ID' for matching
            sequence = seq_d.get(accession_id)
            if sequence:
                row['RNA'] = sequence
            else:
                row['RNA'] = "Sequence_Not_Found"
                not_found_count += 1
            writer.writerow(row)
            if i % 1000 == 0:
                sys.stderr.write(f'Processed {i} rows\n')
    
    sys.stderr.write(f'Total processed rows: {i}\n')
    sys.stderr.write(f'Total sequences not found: {not_found_count}\n')

if __name__ == "__main__":
    main()

