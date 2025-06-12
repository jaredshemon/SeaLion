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
import os
import matplotlib
import matplotlib.pyplot as plt
import re
import subprocess
import numpy as np
from datetime import datetime
from numpy import savetxt
from collections import defaultdict
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import gzip
import matplotlib.patches as mpatches
import pandas as pd

'''
now = datetime.now()
now_format = now.strftime("%Y-%m-%d_%H-%M-%S")


def make_clade_files(fasta_path, clade_path):

    if not os.path.exists(clade_path):
        os.makedirs(clade_path)

    def extract_number(file_name):
        #match = re.search(r'\d+', file_name)
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
                #model_match = re.search(r'Model\s*=\s*(\S+)', line)
                #if model_match:
                 #   model = model_match.group(1)
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
'''
def gc_content(seq):
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    return (g + c) / total * 100 if total > 0 else 0
file_diff = {} 
file_bal = {} 
file_inc = {} 
A_seq_GC = {}
B_seq_GC = {}
C_seq_GC = {}
D_seq_GC = {}
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
                A_seq_GC1 = (round(abs(medians['A']), 3))
                B_seq_GC2 = ((round(abs(medians['B']), 3)))
                C_seq_GC3 = ((round(abs(medians['C']), 3)))
                D_seq_GC4 = ((round(abs(medians['D']), 3)))
                GCinc = (round(abs(medians['B'] + medians['C']), 3))/2
                GCbal = (round(abs(medians['A'] + medians['D']), 3))/2
                GCdif = round(abs(GCinc - GCbal), 3)
                file_diff[file] = GCdif
                file_bal[file] = GCbal
                file_inc[file] = GCinc
                A_seq_GC[file] = A_seq_GC1
                B_seq_GC[file] = B_seq_GC2
                C_seq_GC[file] = C_seq_GC3
                D_seq_GC[file] = D_seq_GC4
        file_diff_sorted_file = dict(sorted(file_diff.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_diff_sorted_GC = dict(sorted(file_diff.items(), key = lambda item: item[1]))
        file_bal_sorted_file = dict(sorted(file_bal.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_inc_sorted_file = dict(sorted(file_inc.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_A_sorted_file = dict(sorted(A_seq_GC.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_B_sorted_file = dict(sorted(B_seq_GC.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_C_sorted_file = dict(sorted(C_seq_GC.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))
        file_D_sorted_file = dict(sorted(D_seq_GC.items(), key=lambda x: int(re.search(r'fastaout(\d+)', x[0]).group(1))))


with open('GC_diff.tsv', 'w') as f:
    f.write(f"GC\tDiff\n")
    for k, v in file_diff_sorted_file.items():
        f.write(f"{k}\t{v}\n")

A_seq_GC_list = []
for k in file_A_sorted_file.values():
    A_seq_GC_list.append(k)

B_seq_GC_list = []
for k in file_B_sorted_file.values():
    B_seq_GC_list.append(k)

C_seq_GC_list = []
for k in file_C_sorted_file.values():
    C_seq_GC_list.append(k)

D_seq_GC_list = []
for k in file_D_sorted_file.values():
    D_seq_GC_list.append(k)


x_axis_bal = []
for k in file_bal_sorted_file.values():
    x_axis_bal.append(k)


x_axis_inc = []
for k in file_inc_sorted_file.values():
    x_axis_inc.append(k)

x_axis = []
for k in file_diff_sorted_file.values():
    x_axis.append(k)

y_axis = range(len(x_axis))


plt.figure(figsize=(16, 6))
plt.plot(y_axis, x_axis, marker='o', linestyle='--', color='blue', label='Δ GC v. AT')

for x in range(1, 60, 10):
    plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

for i in range(1, 60, 20):  # every other bin
    plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

gc_labels = [1, 2, 3, 4, 5, 6]
for i, gc in enumerate(gc_labels):
    plt.text(i * 10 + 5, 24, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
plt.xlabel('Dataset #')
plt.ylabel('Δ GC v. AT')
plt.title('Δ GC v. AT per Dataset')
plt.xticks(ticks=range(len(x_axis)), labels=[i for i in range(len(x_axis))])
plt.yticks(ticks=range(26), labels=[i for i in range(26)])
plt.ylim(0, 26)
plt.legend()
plt.tight_layout()

# Save and show
plt.savefig(f'Δ_GC_v_AT_per_Dataset.svg', dpi=300)
plt.show()                       
plt.close()

plt.figure(figsize=(16, 6))
plt.plot(y_axis, x_axis_inc, marker='o', linestyle='--', color='green', label='GC Increase Clades (B & C)')
plt.plot(y_axis, x_axis_bal, marker='o', linestyle='--', color='orange', label='GC balanced Clades (A & D)')

for x in range(1, 60, 10):
    plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

for i in range(1, 60, 20):  # every other bin
    plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

gc_labels = [1, 2, 3, 4, 5, 6]
for i, gc in enumerate(gc_labels):
    plt.text(i * 10 + 5, 74, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
plt.xlabel('Dataset #')
plt.ylabel('GC Increased and Balanced')
plt.title('GC Content Distribution Across Combined Clades with Bin Highlighting')
plt.ylim(40, 75)
plt.legend()

# Save and show
plt.savefig(f'GC_Increase_v._Balanced_Combined.svg', dpi=300)
plt.show()                       
plt.close()

plt.figure(figsize=(16, 6))
plt.plot(y_axis, A_seq_GC_list, marker='o', linestyle='--', color='blue', label='GC Increase Clade (A)')
plt.plot(y_axis, B_seq_GC_list, marker='o', linestyle='--', color='lightcoral', label='GC Increase Clade (B)')
plt.plot(y_axis, C_seq_GC_list, marker='o', linestyle='--', color='green', label='GC Increase Clade (C)')
plt.plot(y_axis, D_seq_GC_list, marker='o', linestyle='--', color='orange', label='GC balanced Clade (D)')

for x in range(1, 60, 10):
    plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

for i in range(1, 60, 20):  # every other bin
    plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

gc_labels = [1, 2, 3, 4, 5, 6]
for i, gc in enumerate(gc_labels):
    plt.text(i * 10 + 5, 74, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
plt.xlabel('Dataset #')
plt.ylabel('GC Increased and Balanced')
plt.title('GC Content Distribution Across Clades with Bin Highlighting')
plt.ylim(40, 75)
plt.legend()

# Save and show
plt.savefig(f'GC_Increase_v._Balanced.svg', dpi=300)
plt.show()                       
plt.close()
'''
'''
import os
import re
import statistics
from collections import defaultdict
import matplotlib.pyplot as plt

def gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if seq else 0

def get_medians_by_prefix(gc_list, prefixes):
    medians = {}
    for prefix in prefixes:
        values = [
            float(val[1])
            for val in gc_list.values()
            if val[0].startswith(prefix)
        ]
        if values:
            medians[prefix] = statistics.median(values)
    return medians

def extract_number(filename):
    match = re.search(r'fastaout(\d+)', filename)
    return int(match.group(1)) if match else 0

# Containers
file_diff = {}
file_bal = {}
file_inc = {}
clade_gc = {k: {} for k in 'ABCD'}

path = '/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-15_20-59-58_10000'
prefixes = ['A', 'B', 'C', 'D']

for fname in os.listdir(path):
    if fname.endswith('fas'):
        file_path = os.path.join(path, fname)
        with open(file_path, 'r') as f:
            headers, sequences, current_seq = [], [], ''
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        sequences.append(current_seq)
                        current_seq = ''
                    headers.append(line[1:])
                else:
                    current_seq += line
            if current_seq:
                sequences.append(current_seq)

        gc_list = {}
        for header, seq in zip(headers, sequences):
            gc = gc_content(seq)
            gc_list[header] = (header, f"{gc:.2f}", len(seq))

        medians = get_medians_by_prefix(gc_list, prefixes)
        filename_key = re.search(r'fastaout\d+\.fas', fname).group()

        if all(prefix in medians for prefix in 'ABCD'):
            A, B, C, D = [round(medians[p], 3) for p in 'ABCD']
            GCinc = round((B + C) / 2, 3)
            GCbal = round((A + D) / 2, 3)
            GCdiff = round(abs(GCinc - GCbal), 3)

            file_diff[filename_key] = GCdiff
            file_bal[filename_key] = GCbal
            file_inc[filename_key] = GCinc
            clade_gc['A'][filename_key] = A
            clade_gc['B'][filename_key] = B
            clade_gc['C'][filename_key] = C
            clade_gc['D'][filename_key] = D

# Sorting helper
def sort_dict_by_filename(d):
    return dict(sorted(d.items(), key=lambda x: extract_number(x[0])))

# Sorted outputs
file_diff_sorted_file = sort_dict_by_filename(file_diff)
file_bal_sorted_file = sort_dict_by_filename(file_bal)
file_inc_sorted_file = sort_dict_by_filename(file_inc)
sorted_clade_gc = {k: sort_dict_by_filename(v) for k, v in clade_gc.items()}

# GC lists for plotting
gc_lists = {k: list(v.values()) for k, v in sorted_clade_gc.items()}
x_axis = list(file_diff_sorted_file.values())
x_axis_bal = list(file_bal_sorted_file.values())
x_axis_inc = list(file_inc_sorted_file.values())
y_axis = range(len(x_axis))

def plot_line(y, lines, title, ylabel, filename, label_height, ylim=None):
    plt.figure(figsize=(16, 6))
    for data, label, color in lines:
        plt.plot(y, data, marker='o', linestyle='--', label=label, color=color)

    for x in range(1, 60, 10):
        plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

    for i in range(1, 60, 20):
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    for i, gc in enumerate([1, 2, 3, 4, 5, 6]):
        plt.text(i * 10 + 5, int(label_height), f'{gc}', ha='center', va='top', fontsize=9)

    plt.xlabel('Dataset #')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(ticks=range(len(x_axis)), labels=list(range(len(x_axis))))
    if ylim:
        plt.ylim(*ylim)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()
    plt.close()


plot_line(
    y_axis,
    [(x_axis, 'Δ GC v. AT', 'blue')],
    'Δ GC v. AT per Dataset',
    'Δ GC v. AT',
    'Δ_GC_v_AT_per_Dataset.svg',
    label_height=24,
    ylim=(0, 26)
)

plot_line(
    y_axis,
    [(x_axis_inc, 'GC Increase Clades (B & C)', 'green'),
     (x_axis_bal, 'GC balanced Clades (A & D)', 'orange')],
    'GC Content Distribution Across Combined Clades with Bin Highlighting',
    'GC Increased and Balanced',
    'GC_Increase_v._Balanced_Combined.svg',
    label_height=74,
    ylim=(40, 75)
    
)

plot_line(
    y_axis,
    [(gc_lists['A'], 'GC Clade (A)', 'blue'),
     (gc_lists['B'], 'GC Clade (B)', 'lightcoral'),
     (gc_lists['C'], 'GC Clade (C)', 'green'),
     (gc_lists['D'], 'GC Clade (D)', 'orange')],
    'GC Content Distribution Across Clades with Bin Highlighting',
    'GC Increased and Balanced',
    'GC_Increase_v._Balanced.svg',
    label_height=74,
    ylim=(40, 75)
)



def shrink_nj_tree(nj_tree_location, clade_output_path):
THIS DOESNT WORK, NJ DOESN'T PRODUCE THE CORRECT TREES
#############################################################################################################################################
###This function should take the newick string from our tree file, cross check it with the original, then graphs the correct ones ###########
#### vs the incorrect ones. (IQTREE Neighbor Joining Method) Same as graph correct outputs function  ########################################
#############################################################################################################################################   
    
    newick_files = [f for f in os.listdir(nj_tree_location)]
    
    full_path = [os.path.join(nj_tree_location, f) for f in newick_files if f.endswith('bionj')]
    newick_strings_files = {}
    for file in full_path:
        with open(file, 'r') as f:
            for line in f:
                line1 = re.sub(r":-?\d+\.\d+(?:[eE][+-]?\d+)?", "", line)
                newick_strings_files[file] = line1.strip()
    print(newick_strings_files)
    #### This below replaces the fasta file with the corresponding clade file, so we can run the subprocess command below    
    clade_def = [f for f in os.listdir(clade_output_path) if f.startswith('clade_def')]
    
    updated_newick_strings = {}
    for clade_file in clade_def:
        clade_num_match = re.search(r'clade_def_file_(\d+)', clade_file)
        if clade_num_match:
            clade_num = clade_num_match.group(1)
            for fasta_file, newick in newick_strings_files.items():
                fasta_num_match = re.search(r'fastaout_(\d+).fa.bionj', fasta_file)
                if fasta_num_match:
                    fasta_num = fasta_num_match.group(1)
                    # Match clade number with fasta number
                    if fasta_num == clade_num:
                        updated_newick_strings[clade_file] = newick
                        print(f"Matched Clade {clade_file} with Neighbor joining Tree File {fasta_file}")
             
    clade_file_location = clade_output_path
    results = {}
    print(updated_newick_strings)
    for k, v in updated_newick_strings.items(): ##### All of this below only works SOMETIMES when this script is run in the terminal command line, it works well as an HPC cluster job
        full_clade = os.path.join(clade_file_location, k)
        shutil.move(full_clade, reroot_directory)
        command = f'python3 ESofT.py "{v}" {k}'
        print(command)
        esoft_run = subprocess.run(command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        output = esoft_run.stdout.strip()
        reroot_command = f'./reroot.o "{output}" {outgroup} {k}'
        print(reroot_command)
        reroot_run =  subprocess.run(reroot_command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        reroot_output = reroot_run.stdout.strip()
        results[k] = reroot_output
        dst_path = os.path.join(reroot_directory, k)
        shutil.move(dst_path, clade_file_location)

    sorted_results = dict(sorted(results.items(), key=lambda x: int(re.search(r'clade_def_file_(\d+)', x[0]).group(1))))
    print(sorted_results)

nj_tree_location = f'/home/s36jshem_hpc/sealion/runs/setup3/iq_output_2025-06-03_20-55-25'
shrink_nj_tree(nj_tree_location, clade_output_path, correct_newick_string_user_data, tq_dist_path, saving_location, working_directory)

def graph_nj_correctness():
    
    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    newick_strings = []
    for file in sorted(os.listdir(iq_output), key = extract_number):
        if file.endswith('txt'):
            file_path = os.path.join(iq_output, file)
            with open(file_path, 'r') as f:
                newick_string = f.read().strip()
                newick_strings.append(newick_string)
    user_newick = correct_newick_string_user_data
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_file_path = os.path.join(newick_corrected_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_corrected_path, 'newickfile_user.txt' )
    results = []
    
    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            f.write(stripped_newick(user_newick))

        command = f"{tq_dist_path}/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=tq_dist_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))
        
    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    #gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    gc_content_labels = [f"{1 + i * 1}" for i in range(6)]
    correct_counts = [results[i:i + 10].count(0) for i in range(0, 60, 10)]
    incorrect_counts = [results[i:i + 10].count(1) for i in range(0, 60, 10)]
    percent_counts = [(correct / (correct + incorrect)) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]
    # Plotting the results
    x = range(6)
    y = percent_counts
    #coefficients = np.polyfit(x, y, deg=1)
    #polynomial = np.poly1d(coefficients)
    #regression_line = polynomial(x)

    #This saves the x,y data as a csv for future graph overlays
    csv_output = {x_axis:y_axis for x_axis,y_axis in zip(list(x),y)}
    numpy_csv = np.array(list(csv_output.items()))
    csv_path = f'{saving_location}/IQTREE_Success_GC_Content.csv'
    savetxt(csv_path, numpy_csv, delimiter=',') 

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, linestyle = '--', label = '_nolegend_', color = 'green')

    custom_legend = Line2D([0], [0], linestyle='None', marker='o', color='green', label='IQTREE Correct')

    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.title('Tree Success by GC Content (IQ-TREE)')
    plt.xticks(x, gc_content_labels)
    plt.legend(handles=[custom_legend], loc='upper right', fontsize='x-small')
    
    plt.savefig(f'{saving_location}/IQTREE_Success_GC_Content.svg', format='svg')
    plt.show()

nj_tree_location = '/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-30_17-57-54'
reroot_directory = '/home/s36jshem_hpc/sealion/'
clade_output_path = '/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-30_17-57-54'#make sure to run this before sealion
outgroup = 'D'
shrink_nj_tree(nj_tree_location, clade_output_path)

'''
def overlay_jermiin():
    ############################################################################################################################################
    #### This function should overlay all the IQTREE outputs from jermiin's setup- branch lengths .01, .025, .05, to show why we chose .05 #####
    ############################################################################################################################################
    path_01 = '/home/s36jshem_hpc/sealion/runs/jermiin/plots/2025-06-12_23-55-14 /IQTREE_SUCCESS.csv'
    path_025 = '/home/s36jshem_hpc/sealion/runs/jermiin/plots/2025-06-13_00-01-16/IQTREE_SUCCESS.csv'
    path_05 = '/home/s36jshem_hpc/sealion/runs/jermiin/plots/2025-06-13_00-03-12/IQTREE_SUCCESS.csv'
    df1 = pd.read_csv(path_01)
    df2 = pd.read_csv(path_025)
    df3 = pd.read_csv(path_05)

    x1, y1 = df1.iloc[:, 0], df1.iloc[:, 1]
    x2, y2 = df2.iloc[:, 0], df2.iloc[:, 1]
    x3, y3 = df3.iloc[:, 0], df3.iloc[:, 1]
    df = pd.DataFrame({
        'GC Bin': x1,
        'Correct Topologies (%) IQTREE Internal Branch = .01': y1,
        'Correct Topologies (%) IQTREE Internal Branch = .025': y2,
        'Correct Topologies (%) IQTREE Internal Branch = .05': y3
    })
    #print(df)
    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='Internal Branch = .01')
    plt.plot(x2, y2, marker='o', linestyle='--', color='green', label='Internal Branch = .025')
    plt.plot(x3, y3, marker='o', linestyle='--', color='orange', label='Internal Branch = .05')

    plt.xlabel('GC Bin')
    plt.ylabel('% Correct Topologies')
    plt.title('Overlay of Correct Topology Matches by Branch Length')
    plt.xticks(ticks=range(6), labels=[f"{1 + i * 1}" for i in range(6)])
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'/home/s36jshem_hpc/sealion/runs/jermiin/plots/Overlay_jermiin_2.svg', dpi=300)
    plt.show()

overlay_jermiin()