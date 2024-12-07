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
            if record.id not in names:
                names.add(record.id)
                res[record.id] = str(record.seq)
            if i % 100000 == 0:
                sys.stderr.write(f'Parsed {i} sequences from FASTA\n')
    else:
        for i, record in enumerate(sio.parse(path, file_format), 1):
            res[record.id] = str(record.seq)
            if i % 100000 == 0:
                sys.stderr.write(f'Parsed {i} sequences from FASTA\n')
    return res

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python script_name.py metadata.(csv|txt|tsv) fasta_directory/ > output.(csv|txt|tsv)\n")
        sys.exit(1)
        
    metadata_path = sys.argv[1]
    fasta_directory = sys.argv[2]
    
    # Determine delimiter based on metadata file extension
    _, ext = os.path.splitext(metadata_path)
    ext = ext.lower()
    if ext == '.csv':
        delimiter = ','
    elif ext == '.txt' or ext == '.tsv':
        delimiter = '\t'
    else:
        sys.stderr.write("Metadata file must be .csv, .txt, or .tsv\n")
        sys.exit(1)
        
    # Define gene files with full paths
    gene_files = {
        'aa_S': os.path.join(fasta_directory, 'gene_S.translation.fasta'),
        'aa_ORF8': os.path.join(fasta_directory, 'gene_ORF8.translation.fasta'),
        'aa_ORF7a': os.path.join(fasta_directory, 'gene_ORF7a.translation.fasta'),
        'aa_ORF3a': os.path.join(fasta_directory, 'gene_ORF3a.translation.fasta'),
        'aa_N': os.path.join(fasta_directory, 'gene_N.translation.fasta'),
    }
    
    # Create a dictionary of sequences for each gene
    gene_seq_dicts = {}
    for gene_name, fasta_file in gene_files.items():
        if not os.path.exists(fasta_file):
            sys.stderr.write(f"FASTA file for {gene_name} not found: {fasta_file}\n")
            sys.exit(1)
        sys.stderr.write(f"Reading sequences for {gene_name} from {fasta_file}\n")
        gene_seq_dicts[gene_name] = to_dict(path=fasta_file, unique=True)
    
    not_found_counts = {gene_name: 0 for gene_name in gene_files}
    
    with open(metadata_path, newline='', encoding='utf-8', errors='replace') as metadata_f:
        reader = csv.DictReader(metadata_f, delimiter=delimiter)
        if 'seqName' not in reader.fieldnames:
            sys.stderr.write("Error: 'seqName' column not found in metadata file\n")
            sys.exit(1)
        fieldnames = reader.fieldnames + list(gene_files.keys())
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        
        for i, row in enumerate(reader, 1):
            seqname = row.get('seqName', '')
            if not seqname:
                sys.stderr.write(f"Error: Missing 'seqName' in row {i}\n")
                continue
            for gene_name, seq_dict in gene_seq_dicts.items():
                sequence = seq_dict.get(seqname)
                if sequence:
                    row[gene_name] = sequence
                else:
                    row[gene_name] = 'Sequence_Not_Found'
                    not_found_counts[gene_name] += 1
            writer.writerow(row)
            if i % 1000 == 0:
                sys.stderr.write(f'Processed {i} rows\n')
    
    sys.stderr.write(f'Total processed rows: {i}\n')
    for gene_name, count in not_found_counts.items():
        sys.stderr.write(f'Total sequences not found for {gene_name}: {count}\n')

if __name__ == "__main__":
    main()

