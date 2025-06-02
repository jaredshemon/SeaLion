import json
import random
import re
import shutil 
import os
import subprocess
from datetime import datetime
import logging
import sys
import traceback
#import matplotlib.pyplot as plt
import numpy as np
import textwrap
#from matplotlib.lines import Line2D
from numpy import savetxt
from collections import defaultdict
import re
import statistics

'''
now = datetime.now()
now_format = now.strftime("%Y-%m-%d_%H-%M-%S")


def make_clade_files(fasta_path, clade_path):

    if not os.path.exists(clade_path):
        os.makedirs(clade_path)

    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    files = sorted(
        [f for f in os.listdir(fasta_path) if f.startswith('fasta') and f.endswith('.fas')],
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
    small_chunks = 10
    num_chunks = len(sequences[group]) // chunk_size
    gc_contents = [47 + (i // 1) * 4 for i in range(num_chunks)]

    for chunk_idx in range(num_chunks):
        group_chunks = {group: [] for group in groups}
        label_chunks = {group: [] for group in groups}
    
        for group in groups:
            start = chunk_idx * chunk_size
            end = start + chunk_size
            sublist = sequences[group][start:end]

            labels = [f"{group}{j}" for j in range(start + 1, end + 1)]
            paired = list(zip(labels, sublist))
            random.shuffle(paired)

            for k in range(0, chunk_size, small_chunks):
                paired_chunk = paired[k:k + small_chunks]
                group_chunks[group].append(paired_chunk)
                label_chunks[group].append([label for label, seq in paired_chunk])

        for j in range(10):
            output_file = os.path.join(clade_path, f"clade_file_{gc_contents[chunk_idx]}_{chunk_idx+1}_{j+1}.fas")
            with open(output_file, 'w') as f:
                for group in groups:
                    for label, seq in group_chunks[group][j]:
                        f.write(f">{label}\n{seq}\n")

            clade_def_file = os.path.join(clade_path, f"clade_def_file_{gc_contents[chunk_idx]}_{chunk_idx+1}_{j+1}.txt")
            with open(clade_def_file, 'w') as f:
                for group in groups:
                    f.write(f"{group}, {', '.join(label_chunks[group][j])}\n")
    
    return "Sequences extracted and saved! "

fasta_path = f'/home/s36jshem_hpc/sealion/runs/iq_output_2025-04-24_10-59-45'
clade_path = f'//home/s36jshem_hpc/sealion/runs/clade_files_{now_format}'

make_clade_files(fasta_path, clade_path)

def mv_clade_files(clade_file_path, sealion_location):
    This function moves the files from the clade_path to the sealion location, then runs the sealion script
      
    for f in os.listdir(clade_file_path):
        if f.startswith('clade'):
            src_path = os.path.join(clade_file_path, f)
            dst_path = os.path.join(sealion_location, f)
            shutil.move(src_path, dst_path)

    clade_files = []
    clade_def_files = []

    for f in os.listdir(sealion_location):
        if f.startswith('clade_def'):
            clade_def_files.append(f)
        elif f.startswith('clade_file'):
            clade_files.append(f)
    
    clade_files.sort(key=lambda x: int(re.search(r'*(\d+)*', x).group(1)))
    clade_def_files.sort(key=lambda x: int(re.search(r'*(\d+)*', x).group(1)))

    both_files = list(zip(clade_files, clade_def_files))
    print(both_files)

    try:
        #subprocess.run(['./indelible', {temp_path}], cwd = '/home/s36jshem_hpc/sealion/runs' , check=True)            
        os.chdir(sealion_location)
        command = f"apptainer exec SeaLion_container_dir sealion1.pl -i clade_file_67_52.fas  -p clade_def_file_67_52.txt  -o D -M '10000' -l '10000' -prt 1 -s"
        os.system(command)
    except Exception:
        print("Error loading file properly")

    print(f'{f} processed through SeaLion')

clade_file_path = '/home/s36jshem_hpc/sealion/runs/clade_files_2025-04-24_15-36-37'
sealion_location = '/home/s36jshem_hpc/sealion/sealion_script'
mv_clade_files(clade_file_path, sealion_location)


def make_clade_files(fasta_path, clade_path):

    if not os.path.exists(clade_path):
        os.makedirs(clade_path)

    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    files = sorted(
        [f for f in os.listdir(fasta_path) if f.startswith('fasta') and f.endswith('.fas')],
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
    small_chunks = 10
    num_chunks = len(sequences[group]) // chunk_size
    gc_contents = [47 + (i // 1) * 4 for i in range(num_chunks)]

    for chunk_idx in range(num_chunks):
        group_chunks = {group: [] for group in groups}
        label_chunks = {group: [] for group in groups}
    
        for group in groups:
            start = chunk_idx * chunk_size
            end = start + chunk_size
            sublist = sequences[group][start:end]

            labels = [f"{group}{j}" for j in range(start + 1, end + 1)]
            paired = list(zip(labels, sublist))
            random.shuffle(paired)

            for k in range(0, chunk_size, small_chunks):
                paired_chunk = paired[k:k + small_chunks]
                group_chunks[group].append(paired_chunk)
                label_chunks[group].append([label for label, seq in paired_chunk])

        for j in range(10):
            output_file = os.path.join(clade_path, f"clade_file_{gc_contents[chunk_idx]}_{chunk_idx+1}_{j+1}.fas")
            with open(output_file, 'w') as f:
                for group in groups:
                    for label, seq in group_chunks[group][j]:
                        f.write(f">{label}\n{seq}\n")

            clade_def_file = os.path.join(clade_path, f"clade_def_file_{gc_contents[chunk_idx]}_{chunk_idx+1}_{j+1}.txt")
            with open(clade_def_file, 'w') as f:
                for group in groups:
                    f.write(f"{group}, {', '.join(label_chunks[group][j])}\n")
    
    return "Sequences extracted and saved! "

fasta_path = f'/home/s36jshem_hpc/sealion/runs/iq_output_{now_format}'
clade_path = f'//home/s36jshem_hpc/sealion/runs/clade_files_{now_format}'


def mv_clade_files(clade_file_path, sealion_location):
    This function moves the files from the clade_path to the sealion location, then runs the sealion script
      
    clade_def_files = [f for f in os.listdir(clade_file_path) if f.startswith('clade_def_file') and f.endswith('.txt')]
    clade_files = [f for f in os.listdir(clade_file_path) if f.startswith('clade_file') and f.endswith('.fas')]
    
    clade_files.sort(key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    clade_def_files.sort(key=lambda x: int(re.search(r'(\d+)', x).group(1)))

    if len(clade_def_files) != len(clade_files):
        print("Mismatch between clade files and clade definition files. Please check your inputs.")
        return
    
    for clade_file, clade_def_file in zip(clade_def_files, clade_files):
        clade_file_path = os.path.join(clade_path, clade_file)
        clade_def_file_path = os.path.join(clade_path, clade_def_file)

        try:
            shutil.move(clade_def_file, sealion_location)
            shutil.move(clade_file, sealion_location)
            #subprocess.run(['./indelible', {temp_path}], cwd = '/home/s36jshem_hpc/sealion/runs' , check=True)            
            os.chdir(sealion_location)
            command = f"apptainer exec SeaLion_container_dir sealion1.pl -i {clade_file}  -p {clade_def_file}  -o D -M '10000' -l '10000' -prt 1 -s"
            os.system(command)
        except Exception:
            print("Error loading file properly")

    print(f'{i} processed through SeaLion')

mv_clade_files('home/s36jshem_hpc/sealion/runs/Clade_folder/clade_files_{now_format}', '/home/s36jshem_hpc/sealion/sealion_script')


# Uncomment and fix the following block if you need to rename sequence identifiers

    for i in indel_output_files:
        with open(os.path.join(IQTREE_path, i), 'r') as f:
            lines = f.readlines()

        count = 0
        out_file = os.path.join(IQTREE_path, f"unique_{i}")
        with open(out_file, 'w') as f:
            for line in lines:
                if any(x in line for x in ['>A', '>B', '>C', '>D']):
                    count += 1 
                    line = line.strip() + str(count) + "\n"
                    print('File processed with new fasta labels')
                f.write(line)
    
'''
'''
import re

# Path to your user file

def run_AliSIM():
    ###############################################################################################################################################################
    #### This function should format the command for AliSIM, write the newick into a file, then simulate 60 files, then output them into the correct folders ######
    ###############################################################################################################################################################
    path = '/home/s36jshem_hpc/sealion/sealion_files/sea_lion/user_file_ALI.txt'
    keywords = {
        'Format': ('format', 2),
        'SEQ_LENGTH': ('seq_length', 2),
        'out_file': ('out_file', 2),
        'Num_Alignments': ('num_alignments', 2),
        'Newick String': ('newick_string', 3),
    }
    dict1 = {}
    # Read and parse file
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Model'):
                model_match = re.search(r'Model\s*=\s*(\S+)', line)
                if model_match:
                    model = model_match.group(1)
            else:
                for key, (var_name, index) in keywords.items():
                    if line.startswith(key):
                        value = line.split('=', 1)[-1].strip()
                        dict1[var_name] = value
                        break  # Once matched, move to next line

    format = dict1.get('format')
    newick_string = dict1.get('newick_string')
    seq_length = dict1.get('seq_length')
    out_file = dict1.get('out_file')
    num_alignments = dict1.get('num_alignments')
    print(dict1)
    #### Now I need to write the newick into a file, so I can utilize that file for the generation of the alignments
    
    with open('newick.nwk', 'w') as f:
        f.write(newick_string)

    command = f'iqtree2 --alisim {out_file} -m {model} -t newick.nwk --out-format {format} --length {seq_length} --num-alignments {num_alignments}'
    try:
        subprocess.run(command, shell=True)
    except Exception as e:
            print(f"Error loading file properly: {e}")

run_AliSIM()
'''

def gc_content(seq):
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    return (g + c) / total * 100 if total > 0 else 0
file_diff = {} 
path = f'/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-15_20-59-58_10000'
for f in os.listdir(path):
    if f.endswith('fas'):
        # Read and parse the FASTA file
        sequences = []
        headers = []
        gc_list = defaultdict(list)
        current_seq = ""
        file = os.path.join(path, f)
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        sequences.append(current_seq)
                        current_seq = ""
                    headers.append(line[1:])  # Remove '>'
                else:
                    current_seq += line
            if current_seq:  # Add the last sequence
                sequences.append(current_seq)
        values1 = []   
        prefixes = ['A', 'B', 'C', 'D']
        medians = {}        
        gc_summary = {}
        for header, seq in zip(headers, sequences):
            gc = gc_content(seq)
            file1 = re.search(r'fastaout(\d+).fas', str(f)).group()
            gc_list[file1] = [(header, f"{gc:.2f}", (len(seq)))]
            for prefix in prefixes:
                values1 = []
                for key, value in gc_list.items():
                    if value[0][0].startswith(prefix):
                        values1.append(float(value[0][1]))
                if values1:
                    medians[prefix] = statistics.median(values1)
                gc_summary[file1] = medians  
          
        for file, medians in gc_summary.items():
            if 'A' in medians and 'D' in medians:
                GCinc = (round(abs(medians['B'] + medians['C']), 3))/2
                GCbal = (round(abs(medians['A'] + medians['D']), 3))/2
                GCdif = round(abs(GCinc - GCbal), 3)
                file_diff[file] = GCdif
        file_diff_sorted_file = dict(sorted(file_diff.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_diff_sorted_GC = dict(sorted(file_diff.items(), key = lambda item: item[1]))
        #for i in range(0, len(file_diff_sorted_GC), 10):
with open('GC_diff.tsv', 'w') as f:
    f.write(f"GC\tDiff\n")
    for k, v in file_diff_sorted_file.items():
        f.write(f"{k}\t{v}\n")

                                   
