#!usr/bin/env python3 

import collections
import os
# Path to the FASTA file

import os
import collections

# Define GC file path
gc_file_path = f'/home/s36jshem_hpc/sealion/plots/clade_files_2025-05-15_20-59-58/GC_table'
os.makedirs(gc_file_path, exist_ok=True)  # Ensure directory exists

# Open output file for writing GC contents
output_file = os.path.join(gc_file_path, 'GC_contents.tsv')
with open(output_file, 'w') as f:
    f.write("Header\t,GC Content (%)\t,Total Length\n")  # Write CSV header

# Process each file
for i in range(1, 61):
    path = f'/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-15_20-59-58_10000/fastaout{i}.fas'
    print(f"Processing: {path}")
    
    # Read file
    with open(path, 'r') as f:
        lines = f.readlines()

    # Process sequences
    sequences = []
    headers = []
    current_seq = ""

    for line in lines:
        line = line.strip()
        if line.startswith(">"):  # Header line
            if current_seq:  
                sequences.append(current_seq)  # Save previous sequence
                current_seq = ""
            headers.append(line)  # Store new header
        else:
            current_seq += line  # Append sequence lines

    if current_seq:
        sequences.append(current_seq)  # Add last sequence

    # Compute GC content for each sequence
    with open(output_file, 'a') as f:  # Append to the output file
        for header, seq in zip(headers, sequences):
            count = collections.Counter(seq)
            gc_sum = count['G'] + count['C']
            total_sum = len(seq)
            gc_content = (gc_sum / total_sum) * 100 if total_sum > 0 else 0  # Avoid division by zero

            # Write to output file in CSV format
            f.write(f"{os.path.basename(path)}\t,{header}\t,{gc_content:.2f}\t,{total_sum}\n")


            




