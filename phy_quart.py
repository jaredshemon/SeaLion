'''
fasta_path = '/Users/jaredshemonsky/Downloads/INDELibleV1.03/iq_output_2025-03-28_09-58-06'
sequences = []
files = []

def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

for file in os.listdir(fasta_path):
    if file.startswith('GTR') and file.endswith('.fas'):
            files.append(file)
    
files = sorted(files, key = extract_number)

A1_seq = []
for i in files:
    i = os.path.join(fasta_path, i)
    with open(i, 'r') as file:
        lines = file.readlines()
        for index, line in enumerate(lines):
              if line.strip().startswith('>A1'):
                    if index +1 < len(lines):
                        A1_seq.append(lines[index+1].strip())

A1_clade = '/Users/jaredshemonsky/Downloads/INDELibleV1.03/Clade_folder/A1.txt'
clade_path = '/Users/jaredshemonsky/Downloads/INDELibleV1.03/Clade_folder'
if not os.path.exists(clade_path):
      os.makedirs(clade_path)


with open(A1_clade, 'w') as f:
    for index, lines in enumerate(A1_seq, start = 1):
        f.write(f'A{index}\n' + lines + '\n')   
'''

import os
import re
def make_clade_files(fasta_path, clade_path):

    if not os.path.exists(clade_path):
        os.makedirs(clade_path)

    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    files = sorted(
        [f for f in os.listdir(fasta_path) if f.startswith('GTR') and f.endswith('.fas')],
        key=extract_number)

    groups = ["A", "B", "C", "D"]
    sequences = {group: [] for group in groups}  # Dictionary to store sequences for each group

    for file_name in files:
        file_path = os.path.join(fasta_path, file_name)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for index, line in enumerate(lines):
                for group in groups:  # Check all groups                
                    if line.strip().startswith(f'>{group}'):
                        if index + 1 < len(lines):  # Ensure there's a sequence line after
                            sequences[group].append(lines[index + 1].strip())

    chunk_size = 100
    num_chunks = len(sequences[group])/chunk_size
    
    for i in range(int(num_chunks)):
        output_file = os.path.join(clade_path, f"clade_file_{i+1}.fas")
        
        with open(output_file, 'w') as f:
            for group in groups:
                start = i * chunk_size
                end = start + chunk_size
                sublist = sequences[f'{group}'][start:end]
                for index, seq in enumerate(sublist, start=start + 1):
                    f.write(f">{group}{index}\n{seq}\n") 

        clade_definition_file = os.path.join(clade_path, f"clade_def_file_{i+1}.txt")    
        with open(clade_definition_file, 'w') as f:
            for group, seq_list in sequences.items():
                f.write(f'{group}, ')
                for q, seq in enumerate(sublist, start = 1):
                    if q == start + len(sublist):
                        f.write(f'{group}{q} ')    
                    else:
                        f.write(f'{group}{q}, ')
                f.write('\n')
    
    
    return "Sequences extracted and saved! "

fasta_path = '/Users/jaredshemonsky/Downloads/INDELibleV1.03/iq_output_2025-04-08_11-28-11'
clade_path = '/Users/jaredshemonsky/Downloads/INDELibleV1.03/Clade_folder'

make_clade_files(fasta_path, clade_path)

