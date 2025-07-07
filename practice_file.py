
import os
import re
import subprocess
import textwrap
from collections import defaultdict
from datetime import datetime
import numpy as np
import pandas as pd
from numpy import savetxt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, LogNorm, Normalize
import seaborn as sns

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
'''
def overlay_jermiin():
    ############################################################################################################################################
    #### This function should overlay all the IQTREE outputs from jermiin's setup- branch lengths .01, .025, .05, to show why we chose .05 #####
    #### only works in marvin                                                                                                              #####
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
'''

def percent_rejected_heatmaps():
#######################################################################################
### This graphs a heatmap of the percent rejected trees for each dataset in trs1-4 ####   
#######################################################################################

    clade_file_rejected = defaultdict(list)
    percent_rejected = []
    accepted = []
    now_format = ['2025-06-02_15-17-28', '2025-06-14_00-46-23', '2025-06-03_20-55-25', '2025-06-14_00-52-22']
    for i in range(1,5):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/trs{i}'
        for j in range(1,61):
            tsv_location = f'{folder_location}/runs_dir/clade_files_{now_format[i-1]}/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
            match = re.search(r'testresult_clade_file_\d+', tsv_location)
            if match:
                    clade_file_location = match.group()
            for data_files in os.listdir(tsv_location):
                if data_files.startswith('Q6_SeaLion') and data_files.endswith('rejections_nap.tsv'):
                    data_file = data_files
                    median_file = os.path.join(tsv_location, data_file)
                    with open(median_file, 'r') as f:
                        for line in f:
                            line.strip() 
                            if 'RISK_DIST' in line:
                                parts = line.split()
                                if 'rejected' in parts:
                                    rejected1 = float(parts[4]) / 10000
                                    percent_rejected.append({i:rejected1*100})
                                    clade_file_rejected[clade_file_location] = float(rejected1)
                                elif 'initially':
                                    accepted1 = parts[4]
                                    accepted.append(accepted1)
    print(percent_rejected)
    trs = []  
    for i in range(4):
        start = i * 60
        end = start + 60
        trs.append([
            list(percent_rejected[j].values())[0]
            for j in range(start, end)
        ])
    trs1, trs2, trs3, trs4 = trs 
    print(trs1, trs2, trs3, trs4)
    
    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='coolwarm', annot=False, fmt=".1f", cbar_kws={'label': 'Percent Rejected'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])
    ax.set_title("Percent Rejected for Each Dataset (Trs1–4)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–4)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4'], rotation=0)

    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/percent_rejected_heatmap_trs1-4.svg', dpi=400)
        
#percent_rejected_heatmaps()

def delta_support_heatmaps():
##############################################################################################################
### This should look at the difference between correct topology support before and after filtering heatmap ###   
##############################################################################################################
    rejected1 = []
    differences = []
    now_format = ['2025-06-02_15-17-28', '2025-06-14_00-46-23', '2025-06-03_20-55-25', '2025-06-14_00-52-22']
    for i in range(1,5):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/trs{i}'
        after_support = []
        before_support = []
        rejected = []
        for j in range(1,61):
            tsv_location = f'{folder_location}/runs_dir/clade_files_{now_format[i-1]}/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
            match = re.search(r'testresult_clade_file_\d+', tsv_location)
            if match:
                    clade_file_location = match.group()
            for data_files in os.listdir(tsv_location):
                if data_files.startswith('MQ1_') and data_files.endswith('average_tree_support.tsv'):
                    data_file = data_files
                    median_file = os.path.join(tsv_location, data_file)
                    with open(median_file, 'r') as f:
                        for line in f:
                            try:
                                if 'median' in line and '(D,(C,(A,B)));' in line:
                                    parts = line.split()
                                    after_support.append(parts[6])
                                    before_support.append(parts[3])
                                    rejected.append(0)
                            except IndexError as e:
                                rejected.append(1)
                                after_support.append(0)
                                before_support.append(0)
                                pass
        diffs = [float(after) - float(before) for before, after in zip(before_support, after_support)]
        for diff in diffs:
            differences.append({i:diff}) 
        rejected1.append({i:rejected})
    trs = []  
    for i in range(4):
        start = i * 60
        end = start + 60
        trs.append([
            list(differences[j].values())[0]
            for j in range(start, end)
        ])
    trs1, trs2, trs3, trs4 = trs 
    
    rejected_mask = np.array([list(d.values())[0] for d in rejected1])  # shape (4, 60)
    masked_trs = np.where(rejected_mask == 1, np.nan, trs)

    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='coolwarm', annot=False, fmt=".1f", cbar_kws={'label': 'Support Δ'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])
    ax.set_title("Δ Correct Topology After Filtering vs. Before Filtering (Trs1–4)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–4)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4'], rotation=0)

    for i in range(rejected_mask.shape[0]):  # 4 trs
        for j in range(rejected_mask.shape[1]):  # 60 datasets
            if rejected_mask[i, j] == 1:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color='black', linewidth=0))
    rejected_patch = mpatches.Patch(color='black', label='Rejected')
    ax.legend(
        handles=[rejected_patch], loc='upper left', bbox_to_anchor=(1.2, 1), borderaxespad=0, frameon=False, ncol=1)    
    cbar = ax.collections[0].colorbar

            
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/delta_b4_after_support_trs1-4.svg', dpi=400)
    
#delta_support_heatmaps()

def unfiltered_sup_heatmaps():
#######################################################################################################
###This graphs a heatmap of the support values for unfiltered, and collects the values for filtered ###
#######################################################################################################
    now_format = ['2025-06-02_15-17-28', '2025-06-14_00-46-23', '2025-06-03_20-55-25', '2025-06-14_00-52-22']
    topology_supports= []
    unfiltered_topology_supports = []

    for k in range(1,5):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/trs{k}'
        newick_strings = [] #FILTERED NEWICKS FROM SEALION (if fully rejected == 0)
        newick_strings1 = [] #UNFILTERED NEWICKS FROM SEALION
        for j in range(1, 61):
            tsv_location = f'{folder_location}/runs_dir/clade_files_{now_format[k-1]}/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
            best_newick = ''
            best_sup = 0
            best_sup1 = 0
            for data_files in os.listdir(tsv_location):
                if data_files.startswith('MQ1_'):
                    data_file = data_files
                    avg_support_file = os.path.join(tsv_location, data_file)
                    with open(avg_support_file, 'r') as f:
                        lines = f.readlines()
                        for i in lines:
                            line = i.strip()
                            if 'median' in i:
                                list_lines = line.split('\t')
                                try:
                                    sup = float(list_lines[6])
                                    newick = list_lines[1]
                                    if best_sup is None or sup > best_sup:
                                        best_sup = sup
                                        best_newick = newick
                                except (IndexError, ValueError):
                                    pass
                                # Unfiltered (should always be present)
                                try:
                                    sup1 = float(list_lines[3])
                                    newick1 = list_lines[1]
                                    if best_sup1 is None or sup1 > best_sup1:
                                        best_sup1 = sup1
                                        best_newick1 = newick1
                                except (IndexError, ValueError):
                                    pass
        

            topology_supports.append({k:best_sup})
            if best_newick1:
                newick_strings1.append(best_newick1)
            if best_newick:
                newick_strings.append(best_newick)
            else:
                newick_strings.append(0)

            unfiltered_topology_supports.append({k:best_sup1})

    trs = []  
    for i in range(4):
        start = i * 60
        end = start + 60
        trs.append([
            list(unfiltered_topology_supports[j].values())[0]
            for j in range(start, end)
        ])
    trs1, trs2, trs3, trs4 = trs 

    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='coolwarm', annot=False, fmt=".1f", cbar_kws={'label': 'Unfiltered Support Values'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])
    ax.set_title("Unfiltered Support Values for Each Dataset (Trs1–4)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–4)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4'], rotation=0)
    print("minimum", np.min(trs) )
    
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/unfiltered_heatmap_trs1-4.svg', dpi=400)
    
    return topology_supports, unfiltered_topology_supports

#unfiltered_topology_supports, topology_supports = unfiltered_sup_heatmaps()

def filtered_sup_heatmaps(topology_supports):
    print(topology_supports)
    trs = []  
    for i in range(4):
        start = i * 60
        end = start + 60
        trs.append([
            list(topology_supports[j].values())[0]
            for j in range(start, end)
        ])
    
    trs1, trs2, trs3, trs4 = trs 
    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='coolwarm', annot=False, fmt=".1f", cbar_kws={'label': 'Filtered Support Values'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])
    ax.set_title("Filtered Support Values for Each Dataset (Trs1–4)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–4)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4'], rotation=0)
    print("minimum", np.min(trs) )
    rejected_mask = np.full((4, 60), np.nan)
    for idx, entry in enumerate(topology_supports):
        for trs, value in entry.items():
            row = trs - 1
            col = idx % 60
            rejected_mask[row, col] = value
    masked_trs = np.where(rejected_mask == 0, np.nan, rejected_mask)

    for i in range(rejected_mask.shape[0]):  # 4 trs
        for j in range(rejected_mask.shape[1]):  # 60 datasets
            if rejected_mask[i, j] == 0:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color='black', linewidth=0))

    rejected_patch = mpatches.Patch(color='black', label='Rejected')
    ax.legend(
        handles=[rejected_patch], loc='upper left', bbox_to_anchor=(1.2, 1), borderaxespad=0, frameon=False, ncol=1)    
    cbar = ax.collections[0].colorbar
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/filtered_heatmap_trs1-4.svg', dpi=400)
    
    return topology_supports

#filtered_sup_heatmaps(topology_supports)

def IQ_support_heatmap():
    ###############################################################################################
    ###This should graph the same topology graph as the two above but for the IQTREE analyis  #####
    ###############################################################################################
    user_newick = "(((A,B),C),D);"
    tq_dist_path = "/home/jshemonsky/sealion/tqdist/tqDist-1.0.2/install/bin/"
    sorted_data1 = []
    for i in range(1,5):
        now_format = ['2025-06-02_15-17-28', '2025-06-14_00-46-23', '2025-06-03_20-55-25', '2025-06-14_00-52-22']
        newick_path = f'/home/jshemonsky/sealion/runs/runs/trs{i}'
        IQ_likeli_loc = f'/home/jshemonsky/sealion/runs/runs/trs{i}/iq_output_{now_format[i-1]}'
        for file in os.listdir(IQ_likeli_loc):
            if file.startswith('fastaout') and file.endswith('ckp.gz'):
                os.chdir(IQ_likeli_loc)
                os.system(f'gunzip {file}')

        newick_scores = []
        diff_newick_scores = {}
        for file in os.listdir(IQ_likeli_loc):
            if file.startswith('fastaout') and file.endswith('ckp'): ###Remember to change the name if you change the output
                file = os.path.join(IQ_likeli_loc, file)
                with open(file, 'r') as f:
                    lines = f.read()
                    log_likelihoods = re.findall(r'\d+:\s+(-\d+\.\d+)', lines)
                    top = sorted(log_likelihoods)[0]
                    diff_newick_scores[file] = [top]
        y_labels = []
        x_labels = []
        for file_name, score in diff_newick_scores.items():
            match = re.search(r'fastaout_(\d+)\.fa\.ckp', file_name) ###Remember to change this name if you change the output
            if match:
                dataset_number = int(match.group(1))
                x_labels.append(dataset_number)
                y_labels.append(float(score[0]))
        '''
        def extract_number(file_name):
            match = re.search(r'\d+', file_name)
            return int(match.group()) if match else -1
        newick_corrected_path = f'/home/jshemonsky/sealion/runs/runs/trs{i}/corrected_IQ_newick_output_{now_format[i-1]}'
        newick_strings = []
        for file in sorted(os.listdir(newick_corrected_path), key = extract_number):
            if file.startswith('corrected') and file.endswith('txt'):
                file_path = os.path.join(newick_corrected_path, file)
                with open(file_path, 'r') as f:
                    newick_string = f.read().strip()
                    newick_strings.append(newick_string)
        
        newick_file_path = os.path.join(newick_corrected_path, 'newickfile1.txt')
        user_newick_path = os.path.join(newick_corrected_path, 'newickfile_user.txt' )
        results = []
        
        for i in newick_strings:
            with open(newick_file_path, 'w') as f:
                f.write(i)
            with open(user_newick_path, 'w') as f:
                f.write(user_newick)

            command = f"{tq_dist_path}/quartet_dist {newick_file_path} {user_newick_path}"
            result = subprocess.run(command, cwd=newick_path, shell=True, capture_output = True, text = True)
            output = result.stdout.strip()
            results.append(int(output))
        '''
        sorted_data = dict(zip(y_labels, x_labels))
        sorted_data = dict(sorted(sorted_data.items(), key=lambda item: item[1]))
        sorted_data1.append(sorted_data)

    trs = []  
    trs1, trs2, trs3, trs4 = [
        [list(sorted_data1[k].keys())[j] for j in range(60)]
        for k in range(4)
    ]
    trs = [trs1, trs2, trs3, trs4]    
    
    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='bwr', annot=False, fmt=".1f", cbar_kws={'label': 'IQTREE Support Values'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])
    ax.set_title("IQTREE Support Values for Each Dataset (Trs1–4)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–4)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4'], rotation=0)
    
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/IQ_supports_heatmap_trs1-4.svg', dpi=400)
    
#IQ_support_heatmap()

def IQ_delta_sup():
    ########################################################################################################################################################
    ### This graphs the difference between the best and second best tree topologys from IQTREE indicated by log-likelihood        ##########################
    ########################################################################################################################################################    
    sorted_data1 = []
    diff_newick_scores = []
    for i in range(1,5):
        now_format = ['2025-06-02_15-17-28', '2025-06-14_00-46-23', '2025-06-03_20-55-25', '2025-06-14_00-52-22']
        newick_path = f'/home/jshemonsky/sealion/runs/runs/trs{i}'
        IQ_likeli_loc = f'/home/jshemonsky/sealion/runs/runs/trs{i}/iq_output_{now_format[i-1]}'
        for file in os.listdir(IQ_likeli_loc):
            if file.startswith('fastaout') and file.endswith('ckp.gz'):
                os.chdir(IQ_likeli_loc)
                os.system(f'gunzip {file}')
    
        for file in os.listdir(IQ_likeli_loc):
            if file.startswith('fastaout') and file.endswith('ckp'): 
                file = os.path.join(IQ_likeli_loc, file)
                with open(file, 'r') as f:
                    lines = f.read()
                    log_likelihoods = re.findall(r'\d+:\s+(-\d+\.\d+)', lines)
                    top_two = sorted(log_likelihoods, reverse = True)[:2]  
                    try:
                        diff = (float(top_two[0]) - float(top_two[1]))*(-1)
                        diff_newick_scores.append({i: diff})
                    except Exception as e:
                        print(f"Only 1 possibility")   
    trs = [diff_newick_scores[i:i + 60] for i in range(0, 240, 60)]
    trs = [[list(d.values())[0] for d in group] for group in trs]
    trs1, trs2, trs3, trs4 = trs 

    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='plasma', norm=LogNorm(), annot=False, fmt=".1f", cbar_kws={'label': 'IQTREE Δ Support Values'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])
    ax.set_title("IQTREE Δ Support Values for Each Dataset (Trs1–4)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–4)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4'], rotation=0)
    
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/IQ_delta_supports_heatmap_trs1-4.svg', dpi=400)
    
#IQ_delta_sup()

def delta_sup_v_tot_sup_heatmap(unfiltered_topology_supports):
###########################################################################################################
### This the total support of the newick tree v. the delta support of the best v. the second best #########   
###########################################################################################################
        
    now_format = ['2025-06-02_15-17-28', '2025-06-14_00-46-23', '2025-06-03_20-55-25', '2025-06-14_00-52-22']
    diff = []
    diff1 = []
    for i in range(1,5):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/trs{i}'
 
        for j in range(1,61):
            tsv_location = f'{folder_location}/runs_dir/clade_files_{now_format[i-1]}/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
            match = re.search(r'testresult_clade_file_\d+', tsv_location)
            if match:
                    clade_file_location = match.group()
            for data_files in os.listdir(tsv_location):
                if data_files.startswith('MQ1'):
                    data_file = data_files
                    median_file = os.path.join(tsv_location, data_file)
                    scores = []
                    scores1 = []
                    with open(median_file, 'r') as f:
                        for line in f:
                            line = line.strip()
                            if 'median' in line:
                                parts = line.split()
                                if len(parts) > 3:
                                    try:
                                        score1 = float(parts[3])
                                        scores1.append(score1)
                                    except Exception:
                                        scores1.append(0)
                                if len(parts) > 6:
                                    try:
                                        score = float(parts[6])
                                        scores.append(score)
                                    except Exception:
                                        scores.append(0)
                    if len(scores) >= 2:
                        top_two = sorted(scores, reverse = True)[0:2] #THE ISSUE WITH THIS IS IT'S ALWAYS POSITIVE, BECAUSE IT SORTS THE LARGEST FIRSTs
                        diff.append({i:top_two[0] - top_two[1]}) 
                        scores = []
                    elif len(scores) < 2:
                        diff.append(0)
                    if scores1:
                        top_two1 = sorted(scores1, reverse = True)[0:2]
                        diff1.append({i:top_two[0] - top_two[1]})
                        scores1 = []

    differences = diff
    differencesU = diff1
    supports = []
    labels = []
    trs1 = []
    trs2 = []
    trs3 = []
    trs4 = []

    for i, x in zip(unfiltered_topology_supports, diff1):
        u_values = list(i.values())[0]
        u_keys = list(i.keys())[0]
        d_values = list(x.values())[0]
        d_keys = list(x.keys())[0]
        if u_keys and d_keys == 1:
            trs1.append({u_values:d_values})    
        elif u_keys and d_keys == 2:
            trs2.append({u_values:d_values})      
        elif u_keys and d_keys == 3:
            trs3.append({u_values:d_values})  
        elif u_keys and d_keys == 4:
            trs4.append({u_values:d_values})  
                    
    # Initialize buckets for each region
    regions = {
        'Strong conflict:': [],
        'Moderate conflict:': [],
        'Weak support:': [],
        'Strong support:': [],
        'Inconsistent signal:': [],
        'Moderate signal:': [],
        'Well supported:': [],
        'Reliable inference:': [],
    }

    datasets = [trs1, trs2, trs3, trs4]
    for i, dataset in enumerate(datasets):
        print(dataset)
        for j, d in enumerate(dataset):
            # Each d is a dict with one key:value pair
            x, y = list(d.items())[0]
            label = f"Setup:{i+1} Dataset #{j+1}"
            if y < 0.3:
                if x < 0.4:
                    regions['Strong conflict:'].append(label)
                elif x < 0.6:
                    regions['Moderate conflict:'].append(label)
                elif x < 0.8:
                    regions['Weak support:'].append(label)
                else:
                    regions['Strong support:'].append(label)
            else:
                if x < 0.4:
                    regions['Inconsistent signal:'].append(label)
                elif x < 0.6:
                    regions['Moderate signal:'].append(label)
                elif x < 0.8:
                    regions['Well supported:'].append(label)
                else:
                    regions['Reliable inference:'].append(label)
    for name, points in regions.items():
        print(f"{name} {len(points)} datasets")
        for p in points:
            print("  ", p)    
            
    df = pd.DataFrame({
        'Region:': list(regions.keys()),
        'Datasets': [', '.join(labels) for labels in regions.values()],
    })
    print(df.to_markdown(index=False))

    
    region_colors = {
        'Strong conflict:': 'red',
        'Moderate conflict:': 'orange',
        'Weak support:': 'yellow',
        'Strong support:': 'lightgreen',
        'Inconsistent signal:': 'violet',
        'Moderate signal:': 'gold',
        'Well supported:': 'lime',
        'Reliable inference:': 'green',
    }

    region_labels = list(region_colors.keys())
    color_list = [region_colors[label] for label in region_labels]

    # List of datasets
    datasets = [trs1, trs2, trs3, trs4]  # Each is a list of dicts

    # Step 1: Build a 2D array of region indices
    region_index_grid = []
    for dataset in datasets:
        row = []
        for d in dataset:
            x = list(d.keys())[0]
            y = list(d.values())[0]
            # Apply your region logic:
            if y < 0.3:
                if x < 0.4:
                    region = 'Strong conflict:'
                elif x < 0.6:
                    region = 'Moderate conflict:'
                elif x < 0.8:
                    region = 'Weak support:'
                else:
                    region = 'Strong support:'
            else:
                if x < 0.4:
                    region = 'Inconsistent signal:'
                elif x < 0.6:
                    region = 'Moderate signal:'
                elif x < 0.8:
                    region = 'Well supported:'
                else:
                    region = 'Reliable inference:'
            row.append(region_labels.index(region))
        region_index_grid.append(row)
             
    region_index_grid = np.array(region_index_grid)

    # Step 2: Plot
    plt.figure(figsize=(18, 8))
    cmap = ListedColormap(color_list)
    ax = sns.heatmap(region_index_grid, cmap=cmap, annot=False, fmt=".1f", cbar_kws={'label': 'Support/Difference Regions by Custom Color'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4'])

    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=region_colors[lab], label=lab) for lab in region_labels]
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.15, 1), loc='upper left', borderaxespad=0., ncol = 1)

    plt.tight_layout()
    plt.show()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/delta_support_v_tot_sup_unfiltered_heatmap_trs1-4.svg', dpi=400)
       
#delta_sup_v_tot_sup_heatmap(unfiltered_topology_supports)