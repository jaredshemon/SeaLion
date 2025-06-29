#!/usr/bin/env python3


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
import matplotlib.pyplot as plt
import numpy as np
import textwrap
from matplotlib.lines import Line2D
from numpy import savetxt
import errno
import multiprocessing
from multiprocessing import pool
import time
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import gzip
import matplotlib.patches as mpatches
import statistics
from collections import defaultdict
import csv
import pandas as pd
import plottable

######################################################################################################
######################################################################################################
#### User Inputs: These are locations where you need to input the depencies for your script ##########
###################################################################################################### 
working_directory = f'/home/jshemonsky/sealion/runs/runs/trs3' #This is specified in the bash script, where you'd like all your files to end up
how_many_files = 60 #This is how many files you're running 
correct_newick_string_user_data = "(((A,B),C),D);" #This is the correct newick string
sealion_container_location = '/share/scientific_bin/singularity/containers/SeaLion_container.sif' #This is where your sealion container is
iq_model = "F81" #This is the model for IQTREE, SeaLion must be specified in the Sealion.pl file
reroot_directory = '/home/jshemonsky/sealion' #This is where the reroot executables/files are located
user_txt_path = "/home/jshemonsky/sealion/sea_lion/user_file_ALI.txt" #This is where you input the txt file that has all you're AliSim inputs
tq_dist_path = "/home/jshemonsky/sealion/tqdist/tqDist-1.0.2/install/bin/" #This is where tq_dist/quartet_dist command is located
perl_script_location = '/home/jshemonsky/sealion/singlularity_sealion/sealion1.022_dev-version.pl' #The exact path where your perl script is located
icebreaker_location = '/home/jshemonsky/sealion/singlularity_sealion/' #Directory where icebreaker is located
######################################################################################################
######################################################################################################
######################################################################################################
def timed_log(func, description, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    logging.info(f'{description} took {elapsed:.2f} seconds')
    return result

if not os.path.exists(working_directory):
    os.makedirs(working_directory)

now_format = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
#This below is for the logging folder
log_folder = f'{working_directory}/log_folder_{now_format}'
if not os.path.exists(log_folder):
   os.makedirs(log_folder)

log_file = os.path.join(log_folder, 'user_input.log')
logging.basicConfig(filename = log_file, level=logging.DEBUG,
        format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%d/%b/%y %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    ####################################################
    ###This function places the error in the log file ##
    ####################################################
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.error('Uncaught exception', exc_info = (exc_type, exc_value, exc_traceback))
    print('An error has occured, please check the log file for more information')
    sys.exit(1)

sys.excepthook = handle_exception
def run_AliSIM(user_txt_path, working_directory, ALI_output_directory):
    # Keyword mapping from config
    keywords = {
        'Format': ('format', 2),
        'SEQ_LENGTH': ('seq_length', 2),
        'Out_File_Prefix': ('out_file_prefix', 2),
        'Number_Alignments_per_GC': ('num_aligns_per_gc', 2),
        'Newick_Template': ('newick_template', 3),
        'Total_Alignments': ('total_aligns', 2),
        'ACTG': ('ACTG', 2),
        'A1C1T1G1': ('A1C1T1G1', 2),
        'Invariant_Positions': ('invariant', 2),
        'outgroup': ('outgroup', 2)
    }

    # Read user input config
    config = {}
    model = None

    with open(user_txt_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Model'):
                model_match = re.search(r'Model\s*=\s*(\S+)', line)
                if model_match:
                    model = model_match.group(1)
            else:
                for key, (var_name, _) in keywords.items():
                    if line.startswith(key):
                        config[var_name] = line.split('=', 1)[-1].strip()
                        break
    print(config)
    logging.info(f"Here is the config: {config}")
    # Extract required values
    format = config.get('format')
    logging.info(f"Here is the Format: {format}")
    seq_length = config.get('seq_length')
    logging.info(f"Here is the Sequence Length: {seq_length}")
    out_file_prefix = config.get('out_file_prefix')
    logging.info(f"Here is the Out File Prefix Length: {out_file_prefix}")
    num_aligns_per_gc = int(config.get('num_aligns_per_gc'))
    logging.info(f"Here is the Number of Alignments per GC Content Length: {num_aligns_per_gc}")
    total_aligns = int(config.get('total_aligns'))
    logging.info(f"Here is the Total Number of Alignments output: {total_aligns}")
    invariant = config.get('invariant')
    logging.info(f"Here is the pinv: {invariant}")
    newick_template = config.get('newick_template')
    logging.info(f"Here is the Newick: {newick_template}")
    outgroup = config.get('outgroup')
    logging.info(f"Here is the outgroup: {outgroup}")

    def parse_gc_vals(gc_str):
        return list(map(int, gc_str.split()))

    def create_range(start, end, step):
        if step == 0 or start == end:
            return [start]
        elif (start < end and step > 0) or (start > end and step < 0):
            return list(range(start, end + (1 if step > 0 else -1), step))
        else:
            raise ValueError(f"Incompatible range parameters: start={start}, end={end}, step={step}")


    # Identify all GC model keys in config (e.g., ACTG, A1C1T1G1, A2C2T2G2, etc.)
    gc_keys = [key for key in config.keys() if (key == 'ACTG' or (key.startswith('A') and key.endswith('G1')))]

    gc_ranges = {}
    for key in gc_keys:
        min_val, max_val, step = parse_gc_vals(config[key])
        gc_ranges[key] = create_range(min_val, max_val, step)

    def calculate_freqs(gc_val):
        gc_content = gc_val / 100
        A = T = (1 - gc_content) / 2
        G = C = gc_content / 2
        return f"{A:.2f}/{C:.2f}/{G:.2f}/{T:.2f}"

    if not os.path.exists(ALI_output_directory):
        os.makedirs(ALI_output_directory)    

    for iteration in range(1, total_aligns + 1):
        newick = newick_template  # start fresh each iteration

        for key, values in gc_ranges.items():
            idx = (iteration - 1) // num_aligns_per_gc
            if idx >= len(values):
                idx = len(values) - 1
            gc_val = values[idx]
            freq_str = calculate_freqs(gc_val)

            # Replace placeholder, e.g., "F{ACTG}" or "F{A1C1T1G1}" in the template
            placeholder = f"F{{{key}}}"
            newick = newick.replace(placeholder, f"F{{{freq_str}}}")

        # Replace invariant placeholder, if any
        if invariant is not None:
            newick = newick.replace("I1", invariant)

        # Now save/write/use the newick as needed
        newick_path = os.path.join(working_directory, f'newick{iteration}.nwk')
        with open(newick_path, "w") as nwk_file:
            nwk_file.write(newick + '\n')   

        command = f'iqtree3 --alisim {out_file_prefix}_{iteration} -m {model} -t {newick_path} --out-format {format} --length {seq_length}'
        try:
            subprocess.run(command, cwd=working_directory, shell=True, check=True)
        except Exception as e:
            print(f"Error running iqtree3 on iteration {iteration}: {e}")

        # Move generated files
        for f in os.listdir(working_directory):
            if f.endswith('.fa') or f.startswith('newick'):
                src_path = os.path.join(working_directory, f)
                dst_path = os.path.join(ALI_output_directory, f)
                shutil.move(src_path, dst_path)

    print(f"Generated {total_aligns} Newick trees and corresponding output files.")

    return outgroup, newick_template

def rename_seq_fasta(ALI_output_directory, IQTREE_path, iq_model):
     #####################################################################################
     #### This function takes the output from AliSIM then runs it through IQTREE #########
     ##################################################################################### 
    
    # Move files from indelible output to IQ-TREE input
    os.makedirs(IQTREE_path, exist_ok=True)
    pattern = re.compile(r'fastaout_\d+\.fa$')
    for f in os.listdir(ALI_output_directory):
        if pattern.match(f):
            src_path = os.path.join(ALI_output_directory, f)
            dst_path = os.path.join(IQTREE_path, f)
            shutil.move(src_path, dst_path)

    # Collect files in IQ-TREE input directory
    indel_output_files = []
    for f in os.listdir(IQTREE_path):
        if pattern.match(f):
            indel_output_files.append(f)

    
    # Run IQ-TREE on each file in the IQ-TREE input directory
    indel_out_files = []
    for f in os.listdir(IQTREE_path):
        if f.startswith('fasta') and f.endswith('.fa'):
            indel_out_files.append(f)
    
    for file_name in indel_out_files:
        try:
            file_path = os.path.join(IQTREE_path, file_name)
            command = f"iqtree3 -s {file_path} -m {iq_model} " #This can run a little funny on occasion, you may have to go back into it and manually rerun it 
            subprocess.run(command, cwd=IQTREE_path, shell=True)
        except Exception as e:
            print(f"Error loading file properly: {e}")

def move_tree_files(iqtree_output_path, tree_output_path):
    #################################################################################################################################
    ####This function moves the tree files generated by IQ-TREE to the specified tree_output directory. #############################
    #################################################################################################################################

    if not os.path.exists(tree_output_path):
        os.makedirs(tree_output_path)   
    
    for f in os.listdir(iqtree_output_path):
        if f.endswith('.treefile'):
            src_path = os.path.join(iqtree_output_path, f)
            dst_path = os.path.join(tree_output_path, f)
            shutil.move(src_path, dst_path)

def make_clade_files(fasta_path, clade_path, sealion_directory):
    #################################################################################################################################
    ### This function makes the clade files and clade def files, each file iterates to 10 then resets to match the order from the ###
    ######### I originally had it count up to 600 (10 sets of 60), but it doesn't align with the tree outputs for comparison ########
    #################################################################################################################################

    if not os.path.exists(clade_path):
        os.makedirs(clade_path)

    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    files = sorted(
        [f for f in os.listdir(fasta_path) if f.startswith('fasta') and f.endswith('.fa')],
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

    chunk_size = 10
    num_chunks = max(len(seq_list) for seq_list in sequences.values()) / chunk_size
    for i in range(int(num_chunks)):
        output_file = os.path.join(clade_path, f"clade_file_{i+1}.fas")
        
        with open(output_file, 'w') as f:
            for group in groups:
                start = i * chunk_size
                end = start + chunk_size
                sublist = sequences[group][start:end]
                for index, seq in enumerate(sublist, start=1):  # Restart index from 1
                    f.write(f">{group}{index}\n{seq}\n")

        clade_definition_file = os.path.join(clade_path, f"clade_def_file_{i+1}.txt")
        with open(clade_definition_file, 'w', newline = '\n') as f:
            for group in groups:
                start = i * chunk_size
                end = start + chunk_size
                sublist = sequences[group][start:end]
                f.write(f'{group},')
                for j in range(1, len(sublist) + 1):  # Also restart here
                    if j == len(sublist):
                        f.write(f'{group}{j}')
                    else:
                        f.write(f'{group}{j},')
                f.write('\n')
    
        
    return "Sequences extracted and saved! "

def shrink_newick(newick_string_location, newick_corrected_path, clade_output_path, reroot_directory, outgroup): 
    ################################################################################################################################################################
    ##### This function should shrink the newick string down without the branch lengths, then it will run it through juliannes python script 'ESOFT' then it #######
    ##### will run it through reroot, then write the corrected newick strings into a new file                   ####################################################
    ################################################################################################################################################################
    newick_files = [f for f in os.listdir(newick_string_location)]
    
    full_path = [os.path.join(newick_string_location, f) for f in newick_files]
    
    newick_strings_files = {}
    for file in full_path:
        with open(file, 'r') as f:
            for line in f:
                line1 = re.sub(r":\d+\.\d+", "", line)
                newick_strings_files[file] = line1.strip()

    #### This below replaces the fasta file with the corresponding clade file, so we can run the subprocess command below    
    clade_def = [f for f in os.listdir(clade_output_path) if f.startswith('clade_def')]

    updated_newick_strings = {}
    for clade_file in clade_def:
        clade_num_match = re.search(r'clade_def_file_(\d+)', clade_file)
        if clade_num_match:
            clade_num = clade_num_match.group(1)
            for fasta_file, newick in newick_strings_files.items():
                fasta_num_match = re.search(r'fastaout_(\d+)', fasta_file)
                if fasta_num_match:
                    fasta_num = fasta_num_match.group(1)
                    # Match clade number with fasta number
                    if fasta_num == clade_num:
                        updated_newick_strings[clade_file] = newick
                        print(f"Matched Clade {clade_file} with Fasta {fasta_file}")
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
    newick_path = newick_corrected_path
    if not os.path.exists(newick_path):
        os.makedirs(newick_path)

    for clade_def, rerooted_newick in sorted_results.items():
        number_of_file = re.search(r'clade_def_file_(\d+)', clade_def)
        number_for_file = number_of_file.group(1)
        output_file = os.path.join(newick_path, f'corrected_newick_{number_for_file}.txt')
        with open(output_file, 'w') as f:
            f.write(rerooted_newick.strip())
 
def save_as_csv(x, y, graph_path, saving_name):
    #This saves the x,y data as a csv for future graph overlays
    try:
        os.makedirs(graph_path, exist_ok=True)
        df = pd.DataFrame({'x': x, 'y': y})
        csv_path = os.path.join(graph_path, saving_name)
        df.to_csv(csv_path, index=False)
        print(f"Saved CSV to {csv_path}")
    except Exception as e:
        print(f"Failed to save CSV: {e}")
      
def graph_correct_outputs(newick_corrected_path, correct_newick_string_user_data, tq_dist_path, graph_path, working_directory):
    #############################################################################################################################################
    ###This function should take the newick string from our tree file, cross check it with the original, then graphs the correct ones ###########
    #### vs the incorrect ones.                                    ##############################################################################
    #############################################################################################################################################
    
    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    newick_strings = []
    for file in sorted(os.listdir(newick_corrected_path), key = extract_number):
        if file.startswith('corrected') and file.endswith('txt'):
            file_path = os.path.join(newick_corrected_path, file)
            with open(file_path, 'r') as f:
                newick_string = f.read().strip()
                newick_strings.append(newick_string)
    user_newick = correct_newick_string_user_data
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = working_directory
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    
    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            f.write(stripped_newick(user_newick))


        command = f"{tq_dist_path}/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd = tq_dist_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))
    
    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    #gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    gc_content_labels = [f"{1 + i * 1}" for i in range(6)]
    correct_counts = [results[i:i + 10].count(0) for i in range(0, 60, 10)]
    incorrect_counts = [results[i:i + 10].count(1) for i in range(0, 60, 10)]
    percent_counts = [(correct / (correct + incorrect)*100) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]

    
    # Plotting the results
    x = range(6)
    y = percent_counts
    coefficients = np.polyfit(x, y, deg=1)
    polynomial = np.poly1d(coefficients)
    regression_line = polynomial(x)

    if not os.path.exists(graph_path):
        os.makedirs(graph_path)

    #This saves the x,y data as a csv for future graph overlays
    save_as_csv(x, y, graph_path, 'IQTREE_SUCCESS.csv')

    #formatted_newick_tree = '\n'.join(textwrap.wrap(user_data['tree'], width=40))
    #formatted_newick_model = '\n'.join(textwrap.wrap(user_data['b1'], width=40))

    #legend_label = f'Newick Branch Lengths:\n{formatted_newick_tree}\n\nNewick Model:\n{formatted_newick_model}'
    #custom_legend = Line2D([0], [0], linestyle = 'None', color='green', label = legend_label)

    plt.plot(x, y, linestyle = '--', label = '_nolegend_', color = 'green')

    #plt.plot(x, regression_line, 'r--', label = 'regression line')

    plt.xlabel('GC Bin')
    plt.ylabel('% Tree Success')
    plt.title('Tree Success by GC Content (IQ-TREE)')
    plt.xticks(x, gc_content_labels)
    plt.ylim(0,101)
    #plt.legend(handles = [custom_legend], loc = 'upper right', fontsize = 'x-small')

    
    plt.savefig(os.path.join(graph_path, f'IQTREE_SUCCESS.svg'), format='svg')
    plt.show()


def process_task(args):
    fas_fn, txt_fn, sealion_runs_dst, perl_script, sealion_container_location, clade_output_path, icebreaker_location = args

    # Create working subdirectory for this task
    runs_dir = os.path.join(sealion_runs_dst, fas_fn)
    os.makedirs(runs_dir, exist_ok=True)

    icebreaker_src = os.path.join(f'{icebreaker_location}', "icebreaker.o")
    icebreaker_dst = os.path.join(runs_dir, "icebreaker.o")
    shutil.copy(icebreaker_src, icebreaker_dst)
    
    # Full paths to input files
    fas_src = os.path.join(clade_output_path, fas_fn)
    txt_src = os.path.join(clade_output_path, txt_fn)
    # Copy the input files to the task folder
    try:
        # Move files
        fas_dst = os.path.join(runs_dir, fas_fn)
        txt_dst = os.path.join(runs_dir, txt_fn)
        shutil.move(fas_src, fas_dst)
        shutil.move(txt_src, txt_dst)
        # Change to task-specific directory
        os.chdir(runs_dir)

        # Build and run the command
        cmd = f"singularity exec {sealion_container_location} {perl_script} -i {fas_fn} -p {txt_fn} -o D -M '1000' -l '10000' -prt 3 -tlrisk 0.5 -s"
        print(f"→ Running in {runs_dir}: {cmd}")
        exit_code = os.system(cmd)

        if exit_code != 0:
            raise RuntimeError(f"Command failed with exit code {exit_code}")

    except (OSError, shutil.Error, RuntimeError) as e:
        print(f"[ERROR] An error occurred: {e}")
        
def run_sea(perl_script_location, sealion_container_location, clade_output_path, sealion_runs_dst):
####################################################################################
#####This script moves the clade files, copies SeaLion, and runs parallel jobs #####
####################################################################################

    perl_script =  perl_script_location
    
    # Create the main run directory
    if not os.path.exists(sealion_runs_dst):
        os.makedirs(sealion_runs_dst)
    else:
        print(f"Directory {sealion_runs_dst} already exists. Reusing it.")

    # Find and sort clade files
    clade_def_files = [f for f in os.listdir(clade_output_path) if f.startswith('clade_def_file') and f.endswith('.txt')]
    clade_files = [f for f in os.listdir(clade_output_path) if f.startswith('clade_file') and f.endswith('.fas')]

    clade_def_files.sort(key=lambda x: tuple(map(int, re.search(r'clade_def_file_(\d+)\.txt', x).groups())))
    clade_files.sort(key=lambda x: tuple(map(int, re.search(r'clade_file_(\d+)\.fas', x).groups())))

    if len(clade_def_files) != len(clade_files):
        print("Mismatch between .txt and .fas files")
        return
    
    # Prepare the list of tasks
    tasks = []
    for fas_fn, txt_fn in zip(clade_files, clade_def_files):
        tasks.append((fas_fn, txt_fn, sealion_runs_dst, perl_script, sealion_container_location, clade_output_path, icebreaker_location))

    # Run in parallel
    with multiprocessing.Pool(processes=how_many_files) as pool:
        pool.map(process_task, tasks)

    print(f"Processed {len(tasks)} file pairs through SeaLion (in parallel!)")

##################################################################
### This next section is focused on graphing the results    ######
##################################################################

def parse_gc_content(newick: str):
    pattern = re.compile(r'\(([A-D])\d:0\.01.*?\)\:[\d\.]+\[\&model=F81\+F\{([ACTG\d]+)\}')
    clade_gc_status = {}
    
    for match in pattern.finditer(newick):
        clade, freqs = match.groups()
        if freqs == "ACTG":
            clade_gc_status[clade] = "GC Increased"
        elif freqs == 'A2C2T2G2':
            clade_gc_status[clade] = "GC Decreased"
        else:
            clade_gc_status[clade] = "Balanced"
    
    return clade_gc_status

def diff_visualizations(clade_output_path, saving_location):
#############################################################################################################################################
###This graphs support values on the y axis and the clade file on the x axis, more specifically the most supported topology. The topology ###
#### is indicated by the color in the legend.                                    ############################################################
#############################################################################################################################################   
    topology_supports = {}
    unfiltered_topology_supports = {}
    newick_strings = [] #FILTERED NEWICKS FROM SEALION (if fully rejected == 0)
    newick_strings1 = [] #UNFILTERED NEWICKS FROM SEALION

    for j in range(1, 61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
        match = re.search(r'clade_files_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', tsv_location)
        match1 = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', tsv_location)
        if match:
            clade_file_location = match.group()
            clade_file_time = match1.group()
        best_newick = ''
        best_sup = 0
        best_sup1 = 0
        
        for data_files in os.listdir(tsv_location):
            if data_files.startswith('MQ1'):
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
       

        topology_supports[f'clade_file_{j}'] = [best_sup, best_newick if best_newick else "N/A"]
        if best_newick1:
            newick_strings1.append(best_newick1)
        if best_newick:
            newick_strings.append(best_newick)
        else:
            newick_strings.append(0)

        unfiltered_topology_supports[f'clade_file_{j}'] = [best_sup1, best_newick1]
    # Build color mapping for unique topologies
    unique_topologies = sorted(set(value[1] for value in topology_supports.values()))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # expand if needed
    color_map = {topo: colors[i % len(colors)] for i, topo in enumerate(unique_topologies)}
    # Prepare plotting
    labels = []
    supports = []
    bar_colors = []
    topologies = []


    for clade_file in sorted(topology_supports.keys(), key=lambda x: int(x.split('_')[-1])):
        support, topology = topology_supports[clade_file]
        labels.append(clade_file.replace('clade_file_', ''))
        supports.append(support)
        bar_colors.append(color_map[topology])
        topologies.append(topology)
    #Make a Table:
    d = {
         'Dataset #': labels,
         'SeaLion Best Filtered Support Value': supports,
         'Best Supported Topologies': topologies
         }
    df = pd.DataFrame(data=d)
    df.set_index('Dataset #', inplace=True)
    df['Best Supported Topologies'] = df['Best Supported Topologies'].replace('N/A', 'Rejected')
    df.to_csv(f'{saving_location}/1_Table_Best_Supported_RISK+DIST_Topology.csv')
    print(df)
    
    # Plot
    plt.figure(figsize=(20, 6))
    plt.bar(labels, supports, color=bar_colors)
    threshold = 0.6
    threshold1 = 0.4
    plt.axhline(y=threshold, color='blue', linestyle=':', linewidth=1, label=f'Good Support')
    plt.axhline(y=threshold1, color='blue', linestyle=':', linewidth=1, label=f'Moderate Conflict')
    plt.text(60, 0.7, 'Good Support', color='black', fontsize=10, va='top')
    plt.text(60, 0.5, 'Moderate Conflict', color='black', fontsize=10, va='top')
    plt.text(60, 0.3, 'High Conflict', color='black', fontsize=10, va='top')

    for idx, (support, topology) in enumerate(zip(supports, topologies)):
        if topology == "N/A":
                plt.scatter(idx, support + 0.02, marker="*", color="red", s=200, label="Rejected" if idx == 0 else "")

    # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    # Add legend
    legend_handles = [plt.Line2D([0], [0], color=color_map[topo], lw=4, label=topo) for topo in color_map if topo != "N/A"]
    asterisk_handle = mlines.Line2D([], [], color='none', marker='*', markersize=14, 
                                markerfacecolor='red', label='Unsupported', linestyle='None')
    legend_handles.append(asterisk_handle)
    plt.legend(handles=legend_handles, loc = 'upper center', bbox_to_anchor=(0.5, -.1), ncol=3)


    plt.xlabel('Dataset')
    plt.ylabel('Support Value (Max 1)')
    plt.title('Best Supported RISK+DIST Filtered SeaLion Topology per Clade File (Colored by Topology)')
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

    if not os.path.exists(saving_location):
        os.makedirs(saving_location)

    plt.savefig(f"{saving_location}/1_Best_Supported_RISK+DIST_Topology.svg", dpi=300)
    plt.close()
    return best_newick, best_sup, saving_location, newick_strings1, clade_file_location, clade_file_time, tsv_location, unfiltered_topology_supports, newick_strings

def unfiltered_quartet_supports(unfiltered_topology_supports, saving_location):
################################################################################
###This should graph the same as above albeit the UNFILTERED supports  #########
################################################################################ 
        # Build color mapping for unique topologies
    
    unique_topologies = sorted(set(value[1] for value in unfiltered_topology_supports.values()))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # expand if needed
    color_map = {topo: colors[i % len(colors)] for i, topo in enumerate(unique_topologies)}

    # Prepare plotting
    labels = []
    supports = []
    bar_colors = []
    topologies = []

    for clade_file in sorted(unfiltered_topology_supports.keys(), key=lambda x: int(x.split('_')[-1])):
        support, topology = unfiltered_topology_supports[clade_file]
        labels.append(clade_file.replace('clade_file_', ''))
        supports.append(support)
        bar_colors.append(color_map[topology])
        topologies.append(topology)
    #Make a table:
    d = {
         'Dataset #': labels,
         'SeaLion Best Unfiltered Support Value': supports,
         'Best Supported Topologies': topologies,
     }
    df = pd.DataFrame(data=d)
    df.set_index('Dataset #', inplace=True)
    df.to_csv(f'{saving_location}/2_Table_Best_Supported_Unfiltered_Topology.csv')
    print(df)

    # Plot
    plt.figure(figsize=(20, 6))
    plt.bar(labels, supports, color=bar_colors)

    threshold = 0.6
    threshold1 = 0.4
    plt.axhline(y=threshold, color='blue', linestyle=':', linewidth=1, label=f'Good Support')
    plt.axhline(y=threshold1, color='blue', linestyle=':', linewidth=1, label=f'Moderate Conflict')
    plt.text(60, 0.7, 'Good Support', color='black', fontsize=10, va='top')
    plt.text(60, 0.5, 'Moderate Conflict', color='black', fontsize=10, va='top')
    plt.text(60, 0.3, 'High Conflict', color='black', fontsize=10, va='top')

    for idx, (support, topology) in enumerate(zip(supports, topologies)):
        if topology == "N/A":
            plt.text(idx, support + 0.02, "✱", ha='center', va='bottom', fontsize=16, color='red', fontweight='bold')

        # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    # Add legend
    legend_handles = [plt.Line2D([0], [0], color=color_map[topo], lw=4, label=topo) for topo in color_map]
    plt.legend(handles=legend_handles, loc = 'upper center', bbox_to_anchor=(0.5, -.1), ncol=3)

    plt.xlabel('Dataset')
    plt.ylabel('Support Value (Max 1)')
    plt.title('Best Supported Unfiltered SeaLion Topology per Clade File (Colored by Topology)')
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

    if not os.path.exists(saving_location):
        os.makedirs(saving_location)

    plt.savefig(f"{saving_location}/2_Best_Supported_Unfiltered_Topology.svg", dpi=300)
    plt.close()

def IQ_quartet_supports(IQ_likeli_loc, newick_corrected_path, user_newick, tq_dist_path, newick_path, saving_location):
    ###############################################################################################
    ###This should graph the same topology graph as the two above but for the IQTREE analyis  #####
    ###############################################################################################
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

    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

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

    sorted_data = dict(zip(y_labels, x_labels))
    sorted_data = dict(sorted(sorted_data.items(), key=lambda item: item[1]))
    #x_labels_sorted, y_values_sorted = zip(*sorted_data)
    x = np.arange(1,61)
    # Assign colors based on results (0 = correct, 1 = incorrect)
    colors = ['orange' if result == 0 else 'blue' for result in results]
    
    #Make a table
    d = {
         'Dataset #': sorted_data.values(),
         'IQ Log-Likelihood Score': sorted_data.keys(),
         'Correct(0) or Incorrect(1) from tqDist': results,
     }
    df = pd.DataFrame(data=d)
    df.set_index('Dataset #', inplace=True)
    df.to_csv(f'{saving_location}/3_Table_IQ_Topology_barchart.csv')
    print(df)
    
    plt.figure(figsize=(16, 6))
    plt.bar(sorted_data.values(), sorted_data.keys(), color=colors, align='center')
    #correct_newick = set(newick_strings[i] for i, result in enumerate(results) if result == 0)
    #incorrect_newick = set(newick_strings[i] for i, result in enumerate(results) if result == 1)
    #correct_legend_label = "".join(correct_newick)
    #incorrect_legend_label = "".join(incorrect_newick)
    correct_topology = next(newick_strings[i] for i, result in enumerate(results) if result == 0)
    try:
        incorrect_topology = next(newick_strings[i] for i, result in enumerate(results) if result == 1)
    except StopIteration as e:
        incorrect_topology = None

    # Build legend
    correct_patch = plt.Line2D([0], [0], color='orange', lw=4, label=f'Correct: {correct_topology}')
    incorrect_patch = plt.Line2D([0], [0], color='blue', lw=4, label=f'Incorrect: {incorrect_topology}')
    plt.legend(handles=[correct_patch, incorrect_patch], loc='upper right')

    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, -79900, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)


    plt.xlabel("Dataset Number")
    plt.ylabel("Log-Likelihood Score")
    plt.ylim(-74000,-80000)
    plt.title("Log-Likelihood Scores by Dataset")
    plt.xticks(range(1,61))
    plt.tight_layout()
    plt.show()

    plt.savefig(f"{saving_location}/3_IQ_Topology_barchart.svg", dpi=300)
    plt.close()
    
    return results

def graph_correct_outputs1(newick_strings1, correct_newick, tq_dist_path, saving_location):
###############################################################################################################################################################
###This function should take the UNFILTERED newick string from our SeaLion output file, cross check it with the original, then graphs the correct ones ########
#### vs the incorrect ones.                                    ################################################################################################
###############################################################################################################################################################
    ##newick_path = location of the corrected newick strings
    
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = working_directory
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    

    for i in newick_strings1: #Newick strings contains the unfiltered score newick strings, and comes from the diff_visualizations() function
        with open(newick_file_path, 'w') as f:
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            f.write(stripped_newick(correct_newick))


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
    percent_counts = [(correct / (correct + incorrect)*100) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]
    
    # Plotting the results
    x = range(6)
    y = percent_counts
    #coefficients = np.polyfit(x, y, deg=1)
    #polynomial = np.poly1d(coefficients)
    #regression_line = polynomial(x)
    d = {
         'GC Bin': x,
         '% of Correct Newick Strings (Unfiltered SeaLion)': percent_counts,
     }
    df = pd.DataFrame(data=d)
    df.set_index('GC Bin', inplace=True)
    df.to_csv(f'{saving_location}/4_Table_Tree_Success_GC_Content_SeaLion.csv')
    print(df)
    
    #formatted_newick_tree = '\n'.join(textwrap.wrap(user_data['tree'], width=40))
    #formatted_newick_model = '\n'.join(textwrap.wrap(user_data['b1'], width=40))

    #legend_label = f'Newick Branch Lengths:\n{formatted_newick_tree}\n\nNewick Model:\n{formatted_newick_model}'
    custom_legend = Line2D([0], [0], linestyle='None', marker='o', color='green', label='SeaLion (Unfiltered)')

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, linestyle='--', marker='o', color='green', label='Data Points')

    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.title('Tree Success by GC Content (SeaLion)')
    plt.xticks(x, gc_content_labels)
    plt.legend(handles=[custom_legend], loc='upper right', fontsize='x-small')
    plt.ylim(0,101)
    
    plt.savefig((f'{saving_location}/4_Tree_Success_GC_Content_SeaLion.svg'))
    csv_path = os.path.join(f'{saving_location}/4_Table_Tree_Success_GC_Content_SeaLion.csv')
    plt.show()

    return csv_path, results

def graph_correct_outputs2(newick_strings, correct_newick, tq_dist_path, saving_location):
###############################################################################################################################################################
###This function should take the FILTERED newick string from our SeaLion output file, cross check it with the original, then graphs the correct ones ##########
#### vs the incorrect ones.                                    ################################################################################################
###############################################################################################################################################################
    ##newick_path = location of the corrected newick strings
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = working_directory
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results_filtered = []
    rejected_focused = []
    for i in newick_strings:
        if i: #Newick strings contains the filtered score newick strings, and comes from the diff_visualizations() function
            with open(newick_file_path, 'w') as f:
                f.write(str(i))
            with open(user_newick_path, 'w') as f:
                f.write(correct_newick)

        
            command = f"{tq_dist_path}quartet_dist {newick_file_path} {user_newick_path}"
            result = subprocess.run(command, cwd=tq_dist_path, shell=True, capture_output = True, text = True)
            output = result.stdout.strip()
            results_filtered.append(int(output))
            rejected_focused.append('unrejected')
        else:
            results_filtered.append(int(1))
            rejected_focused.append('rejected')

    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    #gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    gc_content_labels = [f"{1 + i * 1}" for i in range(6)]
    correct_counts = [results_filtered[i:i + 10].count(0) for i in range(0, 60, 10)]
    incorrect_counts = [results_filtered[i:i + 10].count(1) for i in range(0, 60, 10)]
    percent_counts = [(correct / (correct + incorrect)*100) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]
    # Plotting the results
    x = range(6)
    y = percent_counts
    #coefficients = np.polyfit(x, y, deg=1)
    #polynomial = np.poly1d(coefficients)
    #regression_line = polynomial(x)
    d = {
         'GC Bin': x,
         '% of Correct Newick Strings (Filtered SeaLion)': percent_counts,
     }
    df = pd.DataFrame(data=d)
    df.set_index('GC Bin', inplace=True)
    df.to_csv(f'{saving_location}/5_Table_RISKDIST_Tree_Success_GC_Content_SeaLion.csv')
    print(df)
    
    
    #formatted_newick_tree = '\n'.join(textwrap.wrap(user_data['tree'], width=40))
    #formatted_newick_model = '\n'.join(textwrap.wrap(user_data['b1'], width=40))

    #legend_label = f'Newick Branch Lengths:\n{formatted_newick_tree}\n\nNewick Model:\n{formatted_newick_model}'
    custom_legend = Line2D([0], [0], linestyle='None', marker='o', color='green', label='SeaLion (RISK+DIST)')

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, linestyle='--', marker='o', color='green', label='Data Points')

    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.title('Tree Success by GC Content (SeaLion)')
    plt.xticks(x, gc_content_labels)
    plt.legend(handles=[custom_legend], loc='upper right', fontsize='x-small')
    plt.ylim(0,101)
    plt.savefig((f'{saving_location}/5_RISKDIST_Tree_Success_GC_Content_SeaLion.svg'))
    csv_path1 = os.path.join(f'{saving_location}/5_Table_RISKDIST_Tree_Success_GC_Content_SeaLion.csv')

    plt.show()

    return csv_path1, results_filtered, rejected_focused
  
def graph_correct_outputsIQ(newick_corrected_path, correct_newick_string_user_data, tq_dist_path, saving_location, working_directory):
#############################################################################################################################################
###This function should take the newick string from our tree file, cross check it with the original, then graphs the correct ones ###########
#### vs the incorrect ones. (IQTREE) Same as graph correct outputs function                         #########################################
#############################################################################################################################################   
    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1

    newick_strings = []
    for file in sorted(os.listdir(newick_corrected_path), key = extract_number):
        if file.startswith('corrected') and file.endswith('txt'):
            file_path = os.path.join(newick_corrected_path, file)
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
    percent_counts = [(correct / (correct + incorrect)*100) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]
    # Plotting the results
    x = range(6)
    y = percent_counts
    #coefficients = np.polyfit(x, y, deg=1)
    #polynomial = np.poly1d(coefficients)
    #regression_line = polynomial(x)
    d = {
         'GC Bin': x,
         '% of Correct Newick Strings (IQ TREE)': percent_counts,
     }
    df = pd.DataFrame(data=d)
    df.set_index('GC Bin', inplace=True)
    df.to_csv(f'{saving_location}/6_Table_IQTREE_Success_GC_Content.csv')
    print(df)
    

    #This saves the x,y data as a csv for future graph overlays
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, linestyle = '--', label = '_nolegend_', color = 'green')

    custom_legend = Line2D([0], [0], linestyle='None', marker='o', color='green', label='IQTREE Correct')

    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.title('Tree Success by GC Content (IQ-TREE)')
    plt.xticks(x, gc_content_labels)
    plt.legend(handles=[custom_legend], loc='upper right', fontsize='x-small')
    plt.ylim(0,101)
    plt.savefig(f'{saving_location}/6_IQTREE_Success_GC_Content.svg', format='svg')
    plt.show()
    
def correct_incorrect_rejected_filtered(results_filtered, rejected_focused, saving_location):
##################################################################################################################
#### This function should overlay the correct vs incorrect vs rejected from Unfiltered SeaLion          ##########
##################################################################################################################

    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    #gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    gc_content_labels = [f"{1 + i * 1}" for i in range(6)]
    correct_counts = [results_filtered[i:i + 10].count(0) for i in range(0, 60, 10)]
    correct_percent = [int(i)/10*100 for i in correct_counts]
    incorrect_counts = [results_filtered[i:i + 10].count(1) for i in range(0, 60, 10)]
    percent_counts = [(correct / (correct + incorrect)) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]
    rejected_counts = [0 if i == 'unrejected' else 1 for i in rejected_focused]
    rejected_counts1 = [rejected_counts[i:i+10].count(1) for i in range(0,60,10)]
    rejected_percent= [int(i)/10*100 for i in rejected_counts1]
    true_incorrect = [sum(1 for j in range(i, i+10) if results_filtered[j] == 1 and rejected_focused[j] == 'unrejected') for i in range(0, 60, 10)]
    incorrect_percent = [int(i)/10*100 for i in true_incorrect]

    x = range(6)
    plt.figure(figsize=(10, 6))
    plt.plot(x, correct_percent, linestyle='--', marker='o', color='green', label='Correct Topologies')
    plt.plot(x, incorrect_percent, linestyle='--', marker='o', color='darkgoldenrod', label='Incorrect Topologies')
    plt.plot(x, rejected_percent, linestyle='--', marker='o', color='darkviolet', label='Rejected Topologies')
    
    df = pd.DataFrame({
        'GC_bin': x,
        'Correct Topologies (%)': correct_percent,
        'Incorrect Topologies (%)': incorrect_percent,
        'Rejected Topologies (%)': rejected_percent
    })
    csv_path = os.path.join(saving_location, f"7_Table_correct_incorrect_rejected_SeaLion.csv")
    df.to_csv(csv_path, index=False)
    print(df)
    
    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.ylim(0,101)
    plt.title('Tree Success, Rejection, and Failure (SeaLion)')
    plt.xticks(x, gc_content_labels)
    plt.legend()

    
    plt.savefig((f'{saving_location}/7_correct_incorrect_rejected_SeaLion.svg'))
    plt.show()

def overlay_correct(csv_path1, IQ_csv_location, saving_location):
###################################################################################################################################
#### This function should overlay the correct vs incorrect from Unfiltered SeaLion v the correct v incorrect from IQTREE ##########
###################################################################################################################################

    df1 = pd.read_csv(IQ_csv_location)
    df2 = pd.read_csv(csv_path1)

    x1, y1 = df1.iloc[:, 0], df1.iloc[:, 1]
    x2, y2 = df2.iloc[:, 0], df2.iloc[:, 1]
    df = pd.DataFrame({
        'GC Bin': x1,
        'Correct Topologies (%) IQTREE': y1,
        'Correct Topologies (%) SeaLion': y2
    })
    csv_path = os.path.join(saving_location, f"8_Table_Unfiltered_Overlay_Tree_Success.csv")
    df.to_csv(csv_path, index=False)
    print(df)
    
    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='ML(IQTREE)')
    plt.plot(x2, y2, marker='o', linestyle='--', color='green', label='SeaLion')

    plt.xlabel('GC Bin')
    plt.ylabel('% Correct Topologies')
    plt.title('Overlay of Correct Topology Matches by GC Content (Unfiltered)')
    plt.xticks(ticks=range(6), labels=[f"{1 + i * 1}" for i in range(6)])
    plt.legend()
    plt.ylim(0,101)
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/8_Unfiltered_Overlay_Tree_Success.svg', dpi=300)
    plt.show()

    return x1, y1, x2, y2

def overlay_correct2(csv_path2, IQ_csv_location, saving_location):
###################################################################################################################################
#### This function should overlay the correct vs incorrect from Filtered SeaLion v the correct v incorrect from IQTREE ############
###################################################################################################################################

    df1 = pd.read_csv(IQ_csv_location)
    df2 = pd.read_csv(csv_path2)

    x1, y1 = df1.iloc[:, 0], df1.iloc[:, 1]
    x2, y2 = df2.iloc[:, 0], df2.iloc[:, 1]
    
    df = pd.DataFrame({
        'GC Bin': x1,
        'Correct Topologies (%) IQTREE': y1,
        'Correct Topologies (%) SeaLion': y2
    })
    
    csv_path = os.path.join(saving_location, f"9_Table_Overlay_Tree_Success.csv")
    df.to_csv(csv_path, index=False)
    print(df)
    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='ML(IQTREE)')
    plt.plot(x2, y2, marker='o', linestyle='--', color='green', label='SeaLion')

    plt.xlabel('GC Content Group')
    plt.ylabel('% Correct Topologies')
    plt.title('Overlay of Correct Topology Matches by GC Content (RISK+DIST)')
    plt.xticks(ticks=range(6), labels=[f"{1 + i * 1}" for i in range(6)])
    plt.legend()
    plt.ylim(0,101)
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/9_Overlay_Tree_Success.svg', dpi=300)
    plt.show()

    return x1, y1, x2, y2

def diff_graphs(clade_output_path, saving_location):
#################################################################################################################################################
### This graphs the difference between the best and second best tree topologys from SeaLion ONLY POSITIVE BECAUSE THE DICT IS SORTED ############
#################################################################################################################################################
    diff = []
    diff1 = []
    for j in range(1,61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
                    diff.append(top_two[0] - top_two[1]) 
                    scores = []
                elif len(scores) < 2:
                    diff.append(0)
                if scores1:
                    top_two1 = sorted(scores1, reverse = True)[0:2]
                    diff1.append(top_two1[0] - top_two1[1])
                    scores1 = []

    bin_avg = []   
    for i in range(0, len(diff), 10):
        gc_bins = diff[i:i+10]
        bin_average = sum(gc_bins)/len(gc_bins)
        bin_avg.append(bin_average)            

    differences = diff
    differencesU = diff1
    y_axis = range(1, 61)
    
    df = pd.DataFrame({
        'Dataset #': y_axis,
        'Δ Best Top Score v. 2nd Best (Filtered)': diff,
        'Δ Best Top Score v. 2nd Best (Unfiltered)': diff1
    })
    
    csv_path = os.path.join(saving_location, f"10_Table_SeaLion_best_second_Δ.csv")
    df.to_csv(csv_path, index=False)
    print(df)
    
    plt.figure(figsize=(16, 5))
    plt.plot(y_axis, differences, marker='+', linestyle='--', color='red', label='Sealion Support Δ' )
    incomplete_x = [i+1 for i,v in enumerate(differences) if v == 0]
    incomplete_y = [0] * len(incomplete_x)
    plt.scatter(incomplete_x, incomplete_y, marker='*', color='black', s=180, label='Fully Rejected')
    plt.axhline(y=.3, color='blue', linestyle=':', linewidth=1, label=f'High Conflict Below')


        # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    plt.xlabel('Dataset')
    plt.ylabel('Support Δ')
    plt.title('Sealion Δ of Best Topology v. Second Best Topology (RISK+DIST)')
    #plt.xticks(ticks=range(50), labels=[f"{47 + i * 4}%" for i in range(5)])
    plt.xticks(ticks=range(1, 61))
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/10_SeaLion_best_second_Δ.svg', dpi=300)
    plt.show()

    return differences, differencesU
    
def diff_graphs1(differencesU, saving_location):
############################################################################################
### Same as the graph above but for unfiltered data ONLY POSITIVE ##########################
############################################################################################

    y_axis = range(1, 61)
    plt.figure(figsize=(16, 5))
    plt.plot(y_axis, differencesU, marker='+', linestyle='--', color='red', label='Sealion Support Δ' )
    plt.axhline(y=.3, color='blue', linestyle=':', linewidth=1, label=f'High Conflict Below')
    
    df = pd.DataFrame({
        'Dataset #': y_axis,
        'Δ Best Top Score v. 2nd Best (Unfiltered)': differencesU
    })
    
    csv_path = os.path.join(saving_location, f"11_Table_SeaLion_Unfil_best_second_Δ.csv")
    df.to_csv(csv_path, index=False)
    print(df)
    
    # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='lightcoral', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    plt.xlabel('Dataset')
    plt.ylabel('Support Δ')
    plt.title('Sealion Δ of Best Topology v. Second Best Topology (Unfiltered)')
    #plt.xticks(ticks=range(50), labels=[f"{47 + i * 4}%" for i in range(5)])
    plt.xticks(ticks=range(1, 61))
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/11_SeaLion_Unfil_best_second_Δ.svg', dpi=300)
    plt.show()

def diff_tree_correct_v_incorrect(differencesU, results, clade_output_path, saving_location):
############################################################################################
### Shows the delta when the tree is correct v. incorrect ONLY POSITIVE ####################
############################################################################################
    for j in range(1,61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
        newicks_scores = {}
        match = re.search(r'testresult_clade_file_\d+', tsv_location)
        if match:
                clade_file_location = match.group()
        for data_files in os.listdir(tsv_location):
            if data_files.startswith('MQ1'):
                data_file = data_files
                median_file = os.path.join(tsv_location, data_file)
                with open(median_file, 'r') as f:
                    for line in f:
                        line = line.strip().split()
                        try:
                            if 'median' in line:
                                newicks_scores[line[1]] = line[3]
                                sorted_dict = dict(sorted(newicks_scores.items(), key = lambda item: item[1], reverse = True)) ###This is just coding for fun I guess, this whole block is useless
                                #print(sort) #This is unfiltered data
                        except Exception as e:
                            pass

    plt.figure(figsize=(20, 6))
    x_axis = range(1,61)
    y_axis = [float(i) for i in differencesU]
    for x, res  in zip(x_axis, results):
        if res == 0:
            plt.bar(x, y_axis[x-1], color = 'skyblue', label = 'Correct Topology Δ' )
        else:
            plt.bar(x, y_axis[x-1], color = 'seagreen', label = 'Incorrect Topology Δ')

    df = pd.DataFrame({
        'Dataset #': y_axis,
        'Δ Best Top Score v. 2nd Best (Unfiltered)': differencesU
    })
    
    csv_path = os.path.join(saving_location, f"12_Table_Δ_Correct_Topology_v._Incorrect_Topology.csv")
    df.to_csv(csv_path, index=False)
    print(df)

    # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .95, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
    threshold = .3
    plt.axhline(y=threshold, color='lightcoral', linestyle=':', linewidth=1)
    plt.text(61, 0.15, 'High Conflict', color='black', fontsize=11, va='top')


    plt.xlabel('Dataset Index')
    plt.ylabel('Δ of Support Values')
    plt.title('Δ Correct Topology v. Incorrect Topology (Unfiltered)')
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.ylim(0,1)
    legend_patch = mpatches.Patch(color='skyblue', label='Correct Topology')
    legend_patch1 = (mpatches.Patch(color='seagreen', label='Incorrect Topology'))
    plt.legend(handles=[legend_patch, legend_patch1], loc = 'upper center', bbox_to_anchor=(0.5, -.1), ncol=3)
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/12_Δ_Correct_Topology_v._Incorrect_Topology.svg', dpi=300)
    plt.show()

def diff_tree_correct_v_incorrect_filtered(differences, results, clade_output_path, saving_location):
############################################################################################
### Shows the delta when the tree is correct v. incorrect ONLY POSITIVE ####################
############################################################################################
    for j in range(1,61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
        newicks_scores = {}
        match = re.search(r'testresult_clade_file_\d+', tsv_location)
        if match:
                clade_file_location = match.group()
        for data_files in os.listdir(tsv_location):
            if data_files.startswith('MQ1'):
                data_file = data_files
                median_file = os.path.join(tsv_location, data_file)
                with open(median_file, 'r') as f:
                    for line in f:
                        line = line.strip().split()
                        try:
                            if 'median' in line:
                                newicks_scores[line[1]] = line[3]
                                sorted_dict = dict(sorted(newicks_scores.items(), key = lambda item: item[1], reverse = True)) ###This is just coding for fun I guess, this whole block is useless
                                #print(sort) #This is unfiltered data
                        except Exception as e:
                            pass

    plt.figure(figsize=(20, 6))
    x_axis = range(1,61)
    y_axis = [float(i) for i in differences]
    for x, res  in zip(x_axis, results):
        if res == 0:
            plt.bar(x, y_axis[x-1], color = 'skyblue', label = 'Correct Topology Δ' )
        else:
            plt.bar(x, y_axis[x-1], color = 'seagreen', label = 'Incorrect Topology Δ')

    df = pd.DataFrame({
        'Dataset #': y_axis,
        'Δ Best Top Score v. 2nd Best (RISK + DIST)': differences
    })
    
    csv_path = os.path.join(saving_location, f"25_Table_Δ_Correct_Topology_v._Incorrect_Topology_Filtered.csv")
    df.to_csv(csv_path, index=False)
    print(df)

    # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .95, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
    threshold = .3
    plt.axhline(y=threshold, color='lightcoral', linestyle=':', linewidth=1)
    plt.text(61, 0.15, 'High Conflict', color='black', fontsize=11, va='top')


    plt.xlabel('Dataset Index')
    plt.ylabel('Δ of Support Values')
    plt.title('Δ Correct Topology v. Incorrect Topology (RISK + DIST)')
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.ylim(0,1)
    legend_patch = mpatches.Patch(color='skyblue', label='Correct Topology')
    legend_patch1 = (mpatches.Patch(color='seagreen', label='Incorrect Topology'))
    plt.legend(handles=[legend_patch, legend_patch1], loc = 'upper center', bbox_to_anchor=(0.5, -.1), ncol=3)
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/25_Δ_Correct_Topology_v._Incorrect_Topology.svg', dpi=300)
    plt.show()

def diff_graphs2(IQ_likeli_loc, saving_location, differences):
########################################################################################################################################################
### This graphs the difference between the best and second best tree topologys from IQTREE indicated by log-likelihood        ##########################
########################################################################################################################################################

    for file in os.listdir(IQ_likeli_loc):
        if file.startswith('fastaout') and file.endswith('ckp.gz'):
            os.chdir(IQ_likeli_loc)
            os.system(f'gunzip {file}')

    diff_newick_scores = []
    for file in os.listdir(IQ_likeli_loc):
        if file.startswith('fastaout') and file.endswith('ckp'): 
            file = os.path.join(IQ_likeli_loc, file)
            with open(file, 'r') as f:
                lines = f.read()
                log_likelihoods = re.findall(r'\d+:\s+(-\d+\.\d+)', lines)
                top_two = sorted(log_likelihoods, reverse = True)[:2]  
                try:
                    diff = (float(top_two[0]) - float(top_two[1]))*(-1)
                    diff_newick_scores.append(diff)
                except Exception as e:
                    print(f"Only 1 possibility")   
    
    y_1 = range(len(differences))
    y_2 = range(len(diff_newick_scores))
    x_1 = differences
    x_2 = diff_newick_scores

    # Convert to numpy array for easier processing
    diffs = np.array(diff_newick_scores)
    indices = np.arange(1, len(diffs) + 1)

    df = pd.DataFrame({
        'Dataset #': indices,
        'Δ Best Top Score v. 2nd Best (IQTREE)': diffs
    })
    
    csv_path = os.path.join(saving_location, f"13_Table_IQ_best_second_Δ.csv")
    df.to_csv(csv_path, index=False)
    print(df)
    
    # Plot with log scaling to highlight small differences
    plt.figure(figsize=(16, 6))
    plt.plot(indices, diffs, marker='+', linestyle='--', color='red', label='IQTREE ΔLnl ' )

    # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, 1e-4, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
    # Highlight significant differences
    threshold = 1e-1
    significant = diffs >= threshold
    #plt.axhline(y=threshold, color='lightcoral', linestyle=':', linewidth=1, label=f'Δ ≥ {threshold}')

    plt.xlabel('Dataset')
    plt.ylabel('Δ Log-Likelihood (Best - 2nd Best)')
    plt.title('ΔLnl (IQ-TREE)')
    plt.xticks(indices, rotation=90)
    plt.yscale('log')  # log scale to emphasize small differences
    plt.ylim(bottom=1e-5, top=max(diffs)*1.5)
    plt.legend()
    #plt.tight_layout()
    plt.savefig(f'{saving_location}/13_IQ_best_second_Δ.svg', dpi=300)
    plt.show()

    return diffs, indices

def diff_graphs3(IQ_likeli_loc, saving_location, differences, results_IQ):
########################################################################################################################################################
### This graphs the difference between the best and second best tree topologys from IQTREE indicated by log-likelihood        ##########################
########################################################################################################################################################
    results = results_IQ
    for file in os.listdir(IQ_likeli_loc):
        if file.startswith('fastaout') and file.endswith('ckp.gz'):
            os.chdir(IQ_likeli_loc)
            os.system(f'gunzip {file}')

    diff_newick_scores = []
    for file in os.listdir(IQ_likeli_loc):
        if file.startswith('fastaout') and file.endswith('ckp'): 
            file = os.path.join(IQ_likeli_loc, file)
            with open(file, 'r') as f:
                lines = f.read()
                log_likelihoods = re.findall(r'\d+:\s+(-\d+\.\d+)', lines)
                top_two = sorted(log_likelihoods, reverse = True)[:2]  
                try:
                    diff = (float(top_two[0]) - float(top_two[1]))*(-1)
                    diff_newick_scores.append(diff)
                except Exception as e:
                    print(f"Only 1 possibility")   
    
    y_1 = range(len(differences))
    y_2 = range(len(diff_newick_scores))
    x_1 = differences
    x_2 = diff_newick_scores

    # Convert to numpy array for easier processing
    diffs = np.array(diff_newick_scores)
    indices = np.arange(1, len(diffs) + 1)
    
    # Plot with log scaling to highlight small differences
    plt.figure(figsize=(16, 6))
    plt.bar(indices, diffs, color='orange', label='IQTREE ΔLnl ' )

    for x, rej in enumerate(results):
        if rej == 1:
            plt.bar(x+1, diffs[x], color='blue')
          
    # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, 1e-1, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
    # Highlight significant differences
    threshold = 1e-1
    significant = diffs >= threshold
    #plt.axhline(y=threshold, color='lightcoral', linestyle=':', linewidth=1, label=f'Δ ≥ {threshold}')

    plt.xlabel('Dataset')
    plt.ylabel('Δ Log-Likelihood (Best - 2nd Best)')
    plt.title('ΔLnl (IQ-TREE)')
    plt.xticks(indices, rotation=90)
    plt.yscale('log')  # log scale to emphasize small differences
    plt.ylim(bottom=1e-6, top=max(diffs)*1.5)
    plt.legend()
    #plt.tight_layout()
    plt.savefig(f'{saving_location}/24_Bar_IQ_best_second_Δ.svg', dpi=300)
    plt.show()

    return diffs, indices

def combined_graph(differences, differencesU, saving_location):
##########################################################################################
### Combined graph with two y-axes to compare filtered and unfiltered datasets ###########
##########################################################################################

    y_axis = range(1, 61)  

    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(16, 5))

    # First dataset (filtered data) on the left y-axis
    ax1.plot(y_axis, differences, marker='.', linestyle='--', color='darkgoldenrod', label='Filtered Δ Support')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Support Δ (Filtered)', color='darkgoldenrod')
    ax1.tick_params(axis='y', labelcolor='darkgoldenrod')

    # Highlight fully rejected datasets with black stars
    incomplete_x = [i + 1 for i, v in enumerate(differences) if v == 0]
    incomplete_y = [0] * len(incomplete_x)
    ax1.scatter(incomplete_x, incomplete_y, marker='.', color='darkgoldenrod', s=180, label='Fully Rejected')
    plt.axhline(y=.3, color='mediumvioletred', linestyle=':', linewidth=1, label=f'High Conflict Below')

    # Add dashed lines and shaded regions for filtered data
    for x in range(11, 61, 10):
        ax1.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)
    for i in range(1, 61, 20):  # Shaded regions
        ax1.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    # GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        ax1.text(i * 10 + 5, 0.98, f'{gc}', ha='center', va='top', fontsize=9, transform=ax1.transData)

    # Second dataset (unfiltered data) on the right y-axis
    ax2 = ax1.twinx()
    ax2.plot(y_axis, differencesU, marker='.', linestyle='--', color='seagreen', label='Unfiltered Δ Support')
    ax2.set_ylim(0.0, 1.0)
    ax2.set_ylabel('Support Δ (Unfiltered)', color='seagreen')
    ax2.tick_params(axis='y', labelcolor='seagreen')


    # Add legends for both datasets
    fig.legend(loc='upper right', bbox_to_anchor=(0.95, 0.85), bbox_transform=ax1.transAxes)

    # Title and layout
    plt.title('Comparison of Filtered and Unfiltered Δ Support')
    plt.xticks(ticks=range(1, 61))


    # Save and show
    df = pd.DataFrame({
        'Dataset #': y_axis,
        'DifferUnfiltered': differencesU,
        'Y_axis_Filtered': differences
    })
    csv_path = os.path.join(saving_location, f"14_Table_Combined_SeaLion_Delta_Support.csv")
    df.to_csv(csv_path, index=False)    
    print(df)
    plt.savefig(f'{saving_location}/14_Combined_SeaLion_Delta_Support.svg', dpi=300)
    plt.show()

def combined_graph_bar(differences, differencesU, saving_location):
#################################################################################################
### Combined graph with two y-axes to compare filtered and unfiltered datasets Barchart #########           
#################################################################################################

    x_axis = np.arange(1, 61)  # Shared x-axis for both datasets
    bar_width = 0.3  # Width of each bar
    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(20, 6))

    # First dataset (filtered data) on the left y-axis
    x_axis = x_axis - bar_width / 2
    ax1.bar(x_axis+bar_width/2, differences, width=bar_width, color='darkgoldenrod', label='Filtered Δ Support')
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Support Δ (Filtered)', color='darkgoldenrod')
    ax1.tick_params(axis='y', labelcolor='darkgoldenrod')
    plt.axhline(y=.3, color='green', linestyle=':', linewidth=1, label=f'Good Support')
    plt.text(61, 0.15, 'High Conflict', color='black', fontsize=11, va='top')

    for x in range(1, 60, 10):
        plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

    for i in range(1, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .98, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    # Plot the unfiltered dataset (differencesU) as bars on the right y-axis
    ax2 = ax1.twinx()
    ax2.bar(x_axis-bar_width/2, differencesU, width=bar_width, color='seagreen', label='Unfiltered Δ Support')
    ax2.set_ylabel('Support Δ (Unfiltered)', color='seagreen')
    ax2.tick_params(axis='y', labelcolor='seagreen')

    # Add legends for both datasets
    fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False, bbox_transform=ax1.transAxes)

    # Title and layout
    plt.title('Comparison of Filtered and Unfiltered Δ Support')
    plt.xticks(ticks=range(1, 61))
    plt.tight_layout
    # Save and show
    df = pd.DataFrame({
        'Dataset #': x_axis,
        'Delta Support Unfiltered': differencesU,
        'Delta Support Filtered': differences
    })
    csv_path = os.path.join(saving_location, f"15_Barchart_Combined_SeaLion_Delta_Support.csv")
    df.to_csv(csv_path, index=False)  
    print(df)
    plt.savefig(f'{saving_location}/15_Barchart_Combined_SeaLion_Delta_Support.svg', dpi=300)
    plt.show()

def combined_graph_IQ(differences, diffs, indices, saving_location):
############################################################################################
### Combined graph with two y-axes to compare filtered and IQ_likelihood delta datasets  ###
############################################################################################

    y_axis = range(1, 61)  # Shared x-axis for both datasets
    gc_labels = [1, 2, 3, 4, 5, 6]

    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(16, 6))

    # Plot filtered dataset on the left y-axis
    ax1.plot(y_axis, differences, linestyle='--', color='seagreen', label='Filtered Δ Support')
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Support Δ (Filtered)', color='seagreen')
    ax1.tick_params(axis='y', labelcolor='seagreen')

    # Highlight fully rejected datasets with black stars
    incomplete_x = [i + 1 for i, v in enumerate(differences) if v == 0]
    incomplete_y = [0] * len(incomplete_x)
    ax1.scatter(incomplete_x, incomplete_y, marker='.', color='black', s=180, label='Fully Rejected')

    '''    # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, 1e-4, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    '''

    # Plot log-likelihood differences on the right y-axis
    ax2 = ax1.twinx()
    ax2.plot(indices, diffs, linestyle='--', color='red', label='Log-Likelihood Δ')
    ax2.set_ylabel('Δ Log-Likelihood (Best - 2nd Best)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_yscale('log')  # Log scaling for likelihood differences

    # Highlight significant log-likelihood differences
    threshold = 0.01
    significant = diffs >= threshold
    #ax2.scatter(indices[significant], diffs[significant], color='red', marker='*', s=100, label=f'Δ ≥ {threshold}')

    # Add legends for both datasets
    fig.legend(loc='upper right', bbox_to_anchor=(0.87, 0.99), bbox_transform=ax1.transAxes)

    # Title and layout
    plt.title('Comparison of Filtered Δ Support and Log-Likelihood Gaps (ΔL)')
    plt.xticks(ticks=range(1, 61))
    plt.tight_layout()

    # Save and show
    df = pd.DataFrame({
        'Dataset #': indices,
        'Δ SeaLion Filtered Support': differences,
        'Δ IQTREE Log-Like Support': diffs
    })
    csv_path = os.path.join(saving_location, f"16_Table_Combined_Delta_SeaLion_LogLikelihood.csv")
    df.to_csv(csv_path, index=False)    
    print(df)
    plt.savefig(f'{saving_location}/16_Combined_Delta_SeaLion_LogLikelihood.svg', dpi=300)
    plt.show()

def combined_graph_IQ_Unfil(differencesU, diffs, indices, saving_location):
##############################################################################################
### Combined graph with two y-axes to compare unfiltered and IQ_likelihood delta datasets ####
##############################################################################################

    y_axis = range(1, 61)  # Shared x-axis for both datasets
    gc_labels = [1, 2, 3, 4, 5, 6]

    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(16, 6))

    # Plot filtered dataset on the left y-axis
    ax1.plot(y_axis, differencesU, linestyle='--', color='seagreen', label='Unfiltered Δ Support')
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Support Δ (Unfiltered)', color='seagreen')
    ax1.tick_params(axis='y', labelcolor='seagreen')

    '''    # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, 1e-4, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    '''

    # Plot log-likelihood differences on the right y-axis
    ax2 = ax1.twinx()
    ax2.plot(indices, diffs, linestyle='--', color='red', label='Log-Likelihood Δ')
    ax2.set_ylabel('Δ Log-Likelihood (Best - 2nd Best)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_yscale('log')  # Log scaling for likelihood differences

    # Highlight significant log-likelihood differences
    threshold = 0.01
    significant = diffs >= threshold
    #ax2.scatter(indices[significant], diffs[significant], color='red', marker='*', s=100, label=f'Δ ≥ {threshold}')

    # Add legends for both datasets
    fig.legend(loc='upper right', bbox_to_anchor=(0.87, 0.99), bbox_transform=ax1.transAxes)

    # Title and layout
    plt.title('Comparison of Unfiltered Δ Support and Log-Likelihood Gaps (ΔL)')
    plt.xticks(ticks=range(1, 61))
    plt.tight_layout()

    # Save and show
    df = pd.DataFrame({
        'Dataset #': indices,
        'Δ SeaLion Unfiltered Support': differencesU,
        'Δ Log_Like Support': diffs
    })
    csv_path = os.path.join(saving_location, f"17_Table_Combined_Support_Delta_Unfiltered_LogLikelihood.csv")
    df.to_csv(csv_path, index=False) 
    print(df)   
    plt.savefig(f'{saving_location}/17_Combined_Support_Delta_Unfiltered_LogLikelihood.svg', dpi=300)
    plt.show()

def support_b4_af_filtering(clade_output_path, results_filtered, saving_location):
################################################################################################################
### This should look at the difference between correct topology support before and after filtering Barchart ####        
################################################################################################################
    results = results_filtered #Because it looks at the results from FILTERED ANALYSIS not the unfiltered analysis
    after_support = []
    before_support = []
    rejected = []
    for j in range(1,61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
    difference = [float(after) - float(before) for before, after in zip(before_support, after_support)]
    x_axis = range(1,61)
    y_1 = [float(i) for i in difference]
    colors = ['darkorange' if i > 0 else 'darkorange' for i in y_1] #The red colors are coordinated with whether or not the filtering increased or decreased the support (originally skyblue)
    plt.figure(figsize=(12, 6))

    # Plot before and after support
    plt.bar(x_axis, y_1, color=colors, edgecolor='black', label='Difference before v. after filtering')
    label_added = False  # Flag to add the label only once

    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [1, 2, 3, 4, 5, 6]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .1, f'{gc}', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    for x, rej in zip(x_axis, rejected):
        if rej == 1:
            plt.scatter(x, 0, marker="x", color="black", s=100, label="rejected" if not label_added else "")
            label_added = True
    
    for x, res in zip(x_axis, results):
        if res == 1:
            plt.bar(x, y_1[x-1], color='red', edgecolor='black', label='Failed Topologies ')
    x_handle = mlines.Line2D([], [], color='none', marker='x', markersize=10, 
                        markerfacecolor='red', label='Rejected', linestyle='None')

    # Labels and title
    plt.xlabel('Dataset Index')
    plt.ylabel('Support Value')
    plt.title('Δ Correct Topology After Filtering vs. Before Filtering')
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.ylim(-0.05, 0.12)
    legend_patch = mpatches.Patch(color='gray', label='Difference After v. Before filtering')
    legend_patch1 = (mpatches.Patch(color='red', label='Failed Topologies (RISK + DIST)'))
    plt.legend(handles=[legend_patch, x_handle, legend_patch1])  
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/18_Correct_Topology_b4_After_Filtering_Bar.svg', dpi=300)
    plt.show()

def support_b4_af(clade_output_path, saving_location):
##########################################################################################################
### This should look at the correct topology support before and after filtering LineGraph (DONT LOVE) ####     
##########################################################################################################

    after_support = []
    before_support = []
    for j in range(1,61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
                                before_support.append(parts[3])
                                after_support.append(parts[6])
                        except IndexError as e:
                            after_support.append(0)
                            pass
    x_axis = range(1,61)
    y_1 = [float(i) for i in before_support]
    y_2 = [float(i) for i in after_support]

    plt.figure(figsize=(12, 6))

    # Plot before and after support
    plt.plot(x_axis, y_1, marker='o', linestyle='-', color='red', label='Before Filtering')
    plt.plot(x_axis, y_2, marker='o', linestyle='-', color='blue', label='After Filtering (RISK + DIST)')

    # Labels and title
    plt.xlabel('Dataset Index')
    plt.ylabel('Support Value')
    plt.title('Correct Topology Before Filtering vs. After Filtering')

    plt.xticks(ticks=x_axis, labels=[str(i) for i in x_axis], rotation=45)    
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    df = pd.DataFrame({
        'Dataset #': x_axis,
        'Support Before Filtering': y_1,
        'Support After Filtering': y_2
    })
    csv_path = os.path.join(saving_location, f"19_SeaLion_correct_topology_b4_after_filtering.csv")
    df.to_csv(csv_path, index=False) 
    print(df)
    plt.savefig(f'{saving_location}/19_SeaLion_correct_topology_b4_after_filtering.svg', dpi=300)
    plt.show()
                           
def reject_GC(clade_output_path, saving_location):
#######################################################################
### This graphs the rejected trees as a function of the GC contents ###   
#######################################################################

    clade_file_rejected = defaultdict(list)
    percent_rejected = []
    accepted = []
    for j in range(1,61):
        tsv_location = f'{clade_output_path}/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
                                percent_rejected.append(rejected1*100)
                                clade_file_rejected[clade_file_location] = float(rejected1)
                            elif 'initially':
                                accepted1 = parts[4]
                                accepted.append(accepted1)
   
    average_percent = []
    gc_bins = []
    for i in range(0, len(percent_rejected), 10): #THIS WOULD JUST PRINT THE AVERAGES, BUT INSTEAD WE'LL PRINT THE DATASETS ONE BY ONE
        group = percent_rejected[i:i+10]
        gc_bins.append(group)
        bin_average = sum(group)/len(group)
        average_percent.append(bin_average)

    x_axis = range(1, 61)
    y_axis = percent_rejected
    plt.figure(figsize=(16, 6))
    df = pd.DataFrame({
        'Dataset #': x_axis,
        '% Rejected Topologies': y_axis,
    })
    csv_path = os.path.join(saving_location, f"20_Table_Bar_SeaLion_percent_rejected.csv")
    df.to_csv(csv_path, index=False) 
    print(df)
    plt.bar(x_axis, y_axis, color='seagreen', edgecolor='black', label='Perecent of Topologies Rejected')


    # Add dashed lines
    for x in range(11, 61, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(1, 61, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    # Add labels, title, and legend
    plt.xlabel('Dataset')
    plt.ylabel('% Rejected')
    plt.title('% of SeaLions Rejected Topologies per Dataset')
    plt.ylim(0, 100)
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Save and show
    plt.savefig(f'{saving_location}/20_Bar_SeaLion_percent_rejected.svg', dpi=300)
    plt.show()

def gc_graphs_multi_models(path, saving_location, newick_template):
##############################################################################################################################################
### This graphs the medians of Clade A, B, C, D. Averages the medians to get a combined graph of GC increasing sequences/GC balanced #########   
##############################################################################################################################################

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
        match = re.search(r'fastaout_(\d+)', filename)
        return int(match.group(1)) if match else 0

    def mean_gc(medians, group):
        vals = [medians[taxon] for taxon in group if taxon in medians]
        return round(sum(vals) / len(vals), 3) if vals else None

    # Containers
    file_diff = {}
    file_bal = {}
    file_inc = {}
    file_dec = {}
    clade_gc = {k: {} for k in 'ABCD'}
    
    prefixes = ['A', 'B', 'C', 'D']
    
    # This parses the newick string 
    newick_string = newick_template
    gc_status = parse_gc_content(newick_string)
    gc_increased = [k for k, v in gc_status.items() if v == "GC Increased"]
    gc_balanced = [k for k, v in gc_status.items() if v == "Balanced"]
    gc_decreased = [k for k, v in gc_status.items() if v == "GC Decreased"]

    for fname in os.listdir(path):
        if fname.endswith('fa'):
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
            filename_key = re.search(r'fastaout_\d+\.fa', fname).group()
            if all(prefix in medians for prefix in 'ABCD'):
                A, B, C, D = [round(medians[p], 3) for p in 'ABCD']
                GCbal = mean_gc(medians, gc_balanced)
                GCinc = mean_gc(medians, gc_increased)####These two variables switch if you switch the sequences that are impacted
                GCdec = mean_gc(medians, gc_decreased)
                GCdiff = round((GCinc - GCbal), 3) #if the graphs are messed up you can edit them here

                file_diff[filename_key] = GCdiff
                file_bal[filename_key] = GCbal
                file_inc[filename_key] = GCinc
                file_dec[filename_key] = GCdec
                clade_gc['A'][filename_key] = A
                clade_gc['B'][filename_key] = B
                clade_gc['C'][filename_key] = C
                clade_gc['D'][filename_key] = D
    # Sorting helper
    def sort_dict_by_filename(d):
        return dict(sorted(d.items(), key=lambda x: extract_number(x[0])))

    # Sorted outputs
    file_diff_sorted_file = sort_dict_by_filename(file_diff)
    file_dec_sorted_file = sort_dict_by_filename(file_dec)
    file_bal_sorted_file = sort_dict_by_filename(file_bal)
    file_inc_sorted_file = sort_dict_by_filename(file_inc)
    sorted_clade_gc = {k: sort_dict_by_filename(v) for k, v in clade_gc.items()}

    #################################
    rows = []
    for idx, filename in enumerate(file_diff_sorted_file):
        if idx % 10 == 0:
            bin_label = idx // 10 + 1
            #print(f"\n--- GC Bin {bin_label} ---")

        # Extract data
        diff = file_diff_sorted_file.get(filename, 'N/A')
        bal = file_bal_sorted_file.get(filename, 'N/A')
        inc = file_inc_sorted_file.get(filename, 'N/A')
        dec = file_dec_sorted_file.get(filename, 'N/A')
        a = sorted_clade_gc['A'].get(filename, 'N/A')
        b = sorted_clade_gc['B'].get(filename, 'N/A')
        c = sorted_clade_gc['C'].get(filename, 'N/A')
        d = sorted_clade_gc['D'].get(filename, 'N/A')

        rows.append({
            'GC Bin': bin_label,
            'Dataset File': filename,
            'Δ GC': diff,
            f'GC Increased ({", ".join(gc_increased)}) Clades': inc,
            f'GC Balanced ({", ".join(gc_balanced)}) Clades': bal,
            f'GC Decreased ({", ".join(gc_balanced)}) Clades': dec,
            'Clade A': a,
            'Clade B': b,
            'Clade C': c,
            'Clade D': d,
        })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(saving_location, "GC_content_table.csv")
    df.to_csv(csv_path, index=False)
    print(df)

    #################################
    
    # GC lists for plotting
    gc_lists = {k: list(v.values()) for k, v in sorted_clade_gc.items()}
    x_axis = list(file_diff_sorted_file.values())
    x_axis_inc = list(file_inc_sorted_file.values())
    x_axis_bal = list(file_bal_sorted_file.values())
    x_axis_dec = list(file_dec_sorted_file.values())
    y_axis = range(len(x_axis))
    def plot_line(y, lines, title, ylabel, filename, label_height, ylim=None):
        plt.figure(figsize=(16, 6))
        for data, label, color in lines:
            plt.plot(y, data, marker='o', linestyle='--', label=label, color=color)
            csv_name = os.path.splitext(os.path.basename(filename))[0] + ".csv"
            save_as_csv(y, data, saving_location, csv_name)

        for x in range(0, 59, 10):
            plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

        for i in range(0, 59, 20):
            plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

        for i, gc in enumerate([1, 2, 3, 4, 5, 6]):
            plt.text(i * 10 + 5, int(label_height), f'{gc}', ha='center', va='top', fontsize=9)

        plt.xlabel('Dataset #')
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xticks(ticks=range(len(x_axis)), labels=list(range(1, len(x_axis)+1)))
        if ylim:
            plt.ylim(*ylim)
        plt.legend()
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        plt.show()
        plt.close()

    color_map = {'A': 'blue', 'B': 'lightcoral', 'C': 'green', 'D': 'orange'}

    plot_line(
        y_axis,
        [(gc_lists[clade], f'GC Clade ({clade}) - {gc_status[clade]}', color_map[clade]) for clade in "ABCD"],
        'GC Content Distribution Across Clades with Bin Highlighting',
        'GC Increased and Balanced',
        f'{saving_location}/21_GC_Increase_v._Balanced.svg', 
        label_height=79,
        ylim=(20, 80)
    )
    
    plot_line(
        y_axis,
        [(x_axis_inc, f'GC Increased Clades ({", ".join(gc_increased)})', 'orange'),
        (x_axis_bal, f'GC Balanced Clades ({", ".join(gc_balanced)})', 'green'),
        *((x_axis_dec, f'GC Decreased Clades ({", ".join(gc_decreased)})', 'blue')if gc_decreased else [])],
        'GC Content Distribution Across Combined Clades with Bin Highlighting',
        'GC Increased and Decreased',
        f'{saving_location}/22_GC_Increase_v._Balanced_Combined.svg',
        label_height=79,
        ylim=(20, 80)
    )

    plot_line(
        y_axis,
        [(x_axis, f'Δ GC Content ({", ".join(gc_balanced)} vs {", ".join(gc_increased)})', 'blue')],
        'Δ GC Content in the Median of Clades {}/{} vs {}/{}'.format(*gc_balanced, *gc_increased, *gc_decreased),
        f'Δ GC Content ({", ".join(gc_balanced)} vs {", ".join(gc_increased)})',
        f'{saving_location}/23_Δ_GC_per_Dataset.svg',
        label_height= 24,
        ylim = (-10, 25)
    )
 
 
##################################################
#FILE LOCATIONS/VARIABLE INPUTS:##################
##################################################
def main():
    now_format = '2025-06-03_20-55-25' 
    # Define all input/output paths here
    ALI_output_directory = f"{working_directory}/ALI_output_{now_format}"  
    iqtree_output_path = f"{working_directory}/iq_output_{now_format}"
    newick_treefile_output_path = f'{working_directory}/tree_output_{now_format}'
    clade_output_path = f"{working_directory}/runs_dir/clade_files_{now_format}/sealion_runs"
    sealion_final_directory = f"{working_directory}/sealion_final_output"
    fasta_path = iqtree_output_path
    newick_corrected_path = f"{working_directory}/corrected_IQ_newick_output_{now_format}"
    sealion_runs_dst = f'{clade_output_path}'
    graph_saving_location = f"{working_directory}/plots/"
    if not os.path.exists(graph_saving_location):
        os.makedirs(graph_saving_location)
    IQ_csv_location = f'{graph_saving_location}IQTREE_SUCCESS.csv'
    
    newick_template = '((((A1:0.01,A2:0.01,A3:0.01,A4:0.01,A5:0.01,A6:0.01,A7:0.01,A8:0.01,A9:0.01,A10:0.01):.39[&model=F81+F{ACTG}+I{I1}],(B1:0.01,B2:0.01,B3:0.01,B4:0.01,B5:0.01,B6:0.01,B7:0.01,B8:0.01,B9:0.01,B10:0.01):.39[&model=F81+F{A1C1T1G1}+I{I1}]):0.05,(C1:0.01,C2:0.01,C3:0.01,C4:0.01,C5:0.01,C6:0.01,C7:0.01,C8:0.01,C9:0.01,C10:0.01):.44[&model=F81+F{A1C1T1G1}+I{I1}]):0.05,(D1:0.01,D2:0.01,D3:0.01,D4:0.01,D5:0.01,D6:0.01,D7:0.01,D8:0.01,D9:0.01,D10:0.01):.49[&model=F81+F{ACTG}+I{I1}]);'

    #outgroup, newick_template = timed_log(run_AliSIM, 'ALISIM', user_txt_path, working_directory, ALI_output_directory)
    #timed_log(rename_seq_fasta, 'IQTREE rename', ALI_output_directory, iqtree_output_path, iq_model)
    #move_tree_files(iqtree_output_path, newick_treefile_output_path)
    #make_clade_files(fasta_path, clade_output_path, sealion_final_directory)
    #shrink_newick(newick_treefile_output_path, newick_corrected_path, clade_output_path, reroot_directory, outgroup)
    #graph_correct_outputs(newick_corrected_path, correct_newick_string_user_data, tq_dist_path, graph_saving_location, working_directory)
    gc_graphs_multi_models(iqtree_output_path, graph_saving_location, newick_template) 
    #run_sea(sealion_container_location, clade_output_path, sealion_runs_dst)
    #best_newick, best_sup, saving_location, newick_strings1, clade_file_location, clade_file_time, tsv_location, unfiltered_topology_supports, newick_strings = diff_visualizations(clade_output_path, graph_saving_location)
    #unfiltered_quartet_supports(unfiltered_topology_supports, graph_saving_location)
    #results_IQ = IQ_quartet_supports(iqtree_output_path, newick_corrected_path, correct_newick_string_user_data, tq_dist_path, working_directory, graph_saving_location)
    #csv_path1, results = graph_correct_outputs1(newick_strings1, correct_newick_string_user_data, tq_dist_path, graph_saving_location)
    #csv_path2, results_filtered, rejected_focused = graph_correct_outputs2(newick_strings, correct_newick_string_user_data, tq_dist_path, graph_saving_location)
    #graph_correct_outputsIQ(newick_corrected_path, correct_newick_string_user_data, tq_dist_path, graph_saving_location, working_directory)#same as the graph_correct_outputs function
    #correct_incorrect_rejected_filtered(results_filtered, rejected_focused, graph_saving_location)
    #x1, y1, x2, y2 = overlay_correct(csv_path1, IQ_csv_location, graph_saving_location)
    #x1, y1, x2, y2 = overlay_correct2(csv_path2, IQ_csv_location, graph_saving_location)
    #differences, differencesU = diff_graphs(clade_output_path, graph_saving_location)
    #diff_graphs1(differencesU, graph_saving_location)
    #diff_tree_correct_v_incorrect(differencesU, results, clade_output_path, graph_saving_location)
    #diff_tree_correct_v_incorrect_filtered(differences, results_filtered, clade_output_path, graph_saving_location)
    #diffs, indices = diff_graphs2(iqtree_output_path, graph_saving_location, differences)
    #diffs, indices = diff_graphs3(iqtree_output_path, graph_saving_location, differences, results_IQ)
    #combined_graph(differences, differencesU, graph_saving_location)
    #combined_graph_bar(differences, differencesU, graph_saving_location)
    #combined_graph_IQ(differences, diffs, indices, graph_saving_location)
    #combined_graph_IQ_Unfil(differencesU, diffs, indices, graph_saving_location)
    ##support_b4_af_filtering(clade_output_path, results_filtered, graph_saving_location)
    #support_b4_af(clade_output_path, graph_saving_location)
    #reject_GC(clade_output_path, graph_saving_location)

if __name__ == "__main__":
    main()   

###################################








'''
def IQ_correct(IQ_csv_location, saving_location):
    #################################################################
    ### This function just replots the IQTREE correct quartets  #####                                              
    #################################################################
    data1 = np.loadtxt(IQ_csv_location, delimiter=',')

    x1, y1 = data1[:, 0], data1[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='ML(IQTREE)')

    plt.xlabel('GC Content Group')
    plt.ylabel('% Tree Success')
    plt.title('IQTree Success by GC Content')
    plt.xticks(ticks=range(6), labels=[f"{47 + i * 4}%" for i in range(6)])
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Tree_Success_GC_Content_IQTREE.svg', dpi=300)
    plt.show()

    return x1, y1

#IQ_correct(IQ_csv_location, graph_saving_location)

'''