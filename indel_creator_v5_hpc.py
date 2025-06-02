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


def timed_log(func, description, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    logging.info(f'{description} took {elapsed:.2f} seconds')
    return result


now_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

#This below is for the logging folder
log_folder = f'/home/s36jshem_hpc/sealion/runs/log_folder_{now_time}'
if not os.path.exists(log_folder):
    os.makedirs(log_folder)

log_file = os.path.join(log_folder, 'user_input.log')
logging.basicConfig(filename = log_file, level=logging.DEBUG,
                    format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%d/%b/%y %H:%M:%S')

def handle_exception(exc_type, exc_value, exc_traceback):
    '''This function places the error in the log file'''
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

    def create_range(min_val, max_val, step):
        if step == 0 or min_val == max_val:
            return [min_val]
        else:
            return list(range(min_val, max_val + 1, step))

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

        command = f'iqtree2 --alisim {out_file_prefix}_{iteration} -m {model} -t {newick_path} --out-format {format} --length {seq_length}'
        try:
            subprocess.run(command, cwd=working_directory, shell=True, check=True)
        except Exception as e:
            print(f"Error running iqtree2 on iteration {iteration}: {e}")

        # Move generated files
        for f in os.listdir(working_directory):
            if f.endswith('.fa') or f.startswith('newick'):
                src_path = os.path.join(working_directory, f)
                dst_path = os.path.join(ALI_output_directory, f)
                shutil.move(src_path, dst_path)

    print(f"Generated {total_aligns} Newick trees and corresponding output files.")

    return outgroup

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
            command = f"iqtree2 -s {file_path} -m {iq_model} "
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
        if file.endswith('txt'):
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
            print(stripped_newick(i))
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            print(stripped_newick(user_newick))
            f.write(stripped_newick(user_newick))


        command = f"./quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd = tq_dist_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))
    
    print(results)
    print('HERE WE ARE')
    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    correct_counts = [results[i:i + 10].count(0) for i in range(0, 60, 10)]
    print(correct_counts)
    incorrect_counts = [results[i:i + 10].count(1) for i in range(0, 60, 10)]
    print(incorrect_counts)
    percent_counts = [(correct / (correct + incorrect)) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]

    
    # Plotting the results
    x = range(6)
    y = percent_counts
    coefficients = np.polyfit(x, y, deg=1)
    polynomial = np.poly1d(coefficients)
    regression_line = polynomial(x)

    if not os.path.exists(graph_path):
        os.makedirs(graph_path)

    #This saves the x,y data as a csv for future graph overlays
    csv_output = {x_axis:y_axis for x_axis,y_axis in zip(list(x),y)}
    numpy_csv = np.array(list(csv_output.items()))
    csv_path = os.path.join(graph_path, f'IQTREE_SUCCESS.csv')
    savetxt(csv_path, numpy_csv, delimiter=',') 

    #formatted_newick_tree = '\n'.join(textwrap.wrap(user_data['tree'], width=40))
    #formatted_newick_model = '\n'.join(textwrap.wrap(user_data['b1'], width=40))

    #legend_label = f'Newick Branch Lengths:\n{formatted_newick_tree}\n\nNewick Model:\n{formatted_newick_model}'
    #custom_legend = Line2D([0], [0], linestyle = 'None', color='green', label = legend_label)

    plt.plot(x, y, linestyle = '--', label = '_nolegend_', color = 'green')

    #plt.plot(x, regression_line, 'r--', label = 'regression line')

    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.title('Tree Success by GC Content (IQ-TREE)')
    plt.xticks(x, gc_content_labels)
    #plt.legend(handles = [custom_legend], loc = 'upper right', fontsize = 'x-small')

    
    plt.savefig(os.path.join(graph_path, f'IQTREE_SUCCESS.png'), format='png')
    plt.show()


def process_task(args):
    fas_fn, txt_fn, sealion_runs_dst, perl_script, sealion_container_location, clade_output_path = args

    # Create working subdirectory for this task
    runs_dir = os.path.join(sealion_runs_dst, fas_fn)
    os.makedirs(runs_dir, exist_ok=True)

    # Full paths to input files
    fas_src = os.path.join(clade_output_path, fas_fn)
    txt_src = os.path.join(clade_output_path, txt_fn)

    # Copy the input files to the task folder
    fas_dst = os.path.join(runs_dir, fas_fn)
    txt_dst = os.path.join(runs_dir, txt_fn)
    shutil.move(fas_src, fas_dst)
    shutil.move(txt_src, txt_dst)

    # Change to task-specific directory
    os.chdir(runs_dir)

    # Build and run the command
    cmd = f"apptainer exec {sealion_container_location} {perl_script} -i {fas_fn} -p {txt_fn} -o D -M '1000' -l '10000' -prt 3 -tlrisk 0.5 -s"
    print(f"→ Running in {runs_dir}: {cmd}")
    os.system(cmd)

def run_sea(sealion_container_location, clade_output_path, sealion_runs_dst):
####################################################################################
#####This script moves the clade files, copies SeaLion, and runs parallel jobs #####
####################################################################################

    perl_script = os.path.join(sealion_container_location, "opt/SeaLion/sealion1.pl")
 
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
        tasks.append((fas_fn, txt_fn, sealion_runs_dst, perl_script, sealion_container_location, clade_output_path))

    # Run in parallel
    with multiprocessing.Pool(processes=60) as pool:
        pool.map(process_task, tasks)

    print(f"Processed {len(tasks)} file pairs through SeaLion (in parallel!)")
    
##################################################
#FILE LOCATIONS/VARIABLE INPUTS:##################
##################################################
def main():
    now_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    now_format = now_time

    # Define all input/output paths here
    user_txt_path = "/home/s36jshem_hpc/sealion/sealion_files/sea_lion/user_file_ALI.txt" #This is where you input the txt file that has all you're AliSim inputs
    working_directory = "/home/s36jshem_hpc/sealion/runs" #This is where you would like to place all of your files
    ALI_output_directory = f"{working_directory}/ALI_output_{now_format}"  
    iqtree_output_path = f"{working_directory}/iq_output_{now_format}"
    newick_treefile_output_path = f'{working_directory}/tree_output_{now_format}'
    clade_output_path = f"/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_{now_format}"
    sealion_final_directory = "/home/s36jshem_hpc/sealion/sealion_script/runs_dir"
    fasta_path = iqtree_output_path
    iq_model = "F81"
    correct_newick_string_user_data = "(((A,B),C),D);"
    tq_dist_path = "/home/s36jshem_hpc/local/bin/"
    graph_path = f"/home/s36jshem_hpc/sealion/plots/clade_files_{now_format}"
    newick_corrected_path = f"{working_directory}/corrected_IQ_newick_output_{now_format}"
    reroot_directory = '/home/s36jshem_hpc/sealion'
    sealion_container_location = '/home/s36jshem_hpc/sealion/sealion_script/SeaLion_container_dir'
    sealion_runs_dst = f'{clade_output_path}/sealion_runs'
    
    outgroup = timed_log(run_AliSIM, 'ALISIM', user_txt_path, working_directory, ALI_output_directory)
    timed_log(rename_seq_fasta, 'IQTREE rename', ALI_output_directory, iqtree_output_path, iq_model)
    move_tree_files(iqtree_output_path, newick_treefile_output_path)
    make_clade_files(fasta_path, clade_output_path, sealion_final_directory)
    shrink_newick(newick_treefile_output_path, newick_corrected_path, clade_output_path, reroot_directory, outgroup)
    graph_correct_outputs(newick_corrected_path, correct_newick_string_user_data, tq_dist_path, graph_path, working_directory)
    #run_sea(sealion_container_location, clade_output_path, sealion_runs_dst)

if __name__ == "__main__":
    main()   

###################################
