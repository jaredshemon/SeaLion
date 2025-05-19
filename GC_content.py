#!usr/bin/env python3 

import collections

# Path to the FASTA file
path = '/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-15_20-59-58_10000/fastaout15.fas'

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
for header, seq in zip(headers, sequences):
    count = collections.Counter(seq)
    gc_sum = count['G'] + count['C']
    total_sum = len(seq)
    gc_content = (gc_sum / total_sum) * 100 if total_sum > 0 else 0  # Avoid division by zero

    # Print results
    print(f"Header: {header}")
    print(f"GC Content: {gc_content:.2f}%")
    print(f"Total Length: {total_sum}")
    print("-" * 30)


