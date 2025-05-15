#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 12:37:55 2025

@author: jaredshemonsky
"""

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

# Read user data from a file
def read_user_data(user_file_path):
    ##########################################################################################################################################
    ###### This reads the userfile information into a dictionary, along with the models which need to be handled #############################
    ###### separately so that they don't overwrite eachother     #############################################################################
    ##########################################################################################################################################

   
    user_data = {}
    model_counter = 1
    current_model = {}
    with open(user_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and '=' in line:
                key, value = line.split('=', 1)
                key, value = key.strip(), value.strip()
                if key == "model":
                    if current_model:
                        user_data[f"model{model_counter}"] = current_model
                        model_counter += 1
                        current_model = {}
                    current_model["name"] = value
                elif key in ["submodel", "statefreq", "rates", "insertmodel", "insertrate"]:
                    current_model[key] = value
                elif key in 'IQ_model':
                    iq_model = value
                elif key in 'randomseed':
                    if value == "":
                        user_data[key] = random.randint(1, 100000)
                    else:
                       user_data[key]= value 
                else:
                    user_data[key] = value
        if current_model:
            user_data[f"model{model_counter}"] = current_model

    logging.info(f'User Input: {user_data}')
    logging.info(f'Input Model: {user_data}')

    return user_data, iq_model

# Validate user data against the template
def validate_user_data(user_data, json_template):
    ##########################################################################################################################
    ##### This funciion validates that whatever is in the ancestralprint, output, and model_name parameters is in ############
    ##### the json template                                      #############################################################
    ##########################################################################################################################

    valid_ancestralprints = json_template["ancestralprint"]
    valid_outputs = json_template["output"]
    valid_models = json_template["models"].keys()

    # Check ancestralprint
    if 'ancestralprint' in user_data and user_data['ancestralprint'] not in valid_ancestralprints:
        raise ValueError(f"Invalid ancestralprint value: {user_data['ancestralprint']}. Valid options are: {valid_ancestralprints}")

    # Check output
    if 'output' in user_data and user_data['output'] not in valid_outputs:
        raise ValueError(f"Invalid output value: {user_data['output']}. Valid options are: {valid_outputs}")

    # Check models
    model_counter = 1
    while f'model{model_counter}' in user_data:
        model_name = user_data[f'model{model_counter}']["name"]
        if model_name not in valid_models:
            raise ValueError(f"Invalid model name: {model_name}. Valid options are: {list(valid_models)}")
        model_counter += 1

# Update the JSON configuration template with user data
def update_json_structure(json_data, user_data, iteration):
    ###################################################################################################################################
    #### This function places the dictionary: user_data's information into the json template, the model function at the     ###########
    #### bottom is needed so the user can input minGC, maxGC, and stepGC, and get a steadily increasing GC content according ##########
    #### to the min/max/step every 100 iterations               #######################################################################
    ###################################################################################################################################
    # Update the JSON structure with the user data
    if 'TYPE' in user_data:
        json_data['TYPE'][0] = user_data['TYPE']
    
    if 'ancestralprint' in user_data:
        json_data['ancestralprint'][0] = user_data['ancestralprint']
    
    if 'output' in user_data:
        json_data['output'][0] = user_data['output']
    
    if 'randomseed' in user_data:
        json_data['randomseed'] = user_data['randomseed']
    
    models = []
    model_counter = 1
    model_counts = {}
    while f'model{model_counter}' in user_data:
        model_data = user_data[f'model{model_counter}']
        model_name = model_data["name"]
        if model_name in model_counts:
            model_counts[model_name] += 1
            model_name = f"{model_name}_{model_counts[model_name]}"
        else:
            model_counts[model_name] = 1    
        
        if isinstance(model_data['statefreq'], str):
            model_data['statefreq'] = model_data['statefreq'].split(' ')
            
        gc_min = int(model_data['statefreq'][0])
        gc_max = int(model_data['statefreq'][1])
        gc_step = int(model_data['statefreq'][2])
        
        if gc_step == 0:
         gc_content = gc_min / 100
         
        else:
            #step_increase = (iteration - 1) // 100
            #gc_content = min(gc_min + step_increase * gc_step, gc_max) / 100
            for i in range(gc_min, gc_max + 1, gc_step):
                if (iteration - 1) // 10 == (i - gc_min) // gc_step:
                    gc_content = i / 100
                    break
     
        A = T = (1 - gc_content) / 2
        G = C = gc_content / 2
        statefreq = f"{T:.2f} {C:.2f} {A:.2f} {G:.2f}"
                                        
        models.append({
            "name": model_name,
            "submodel": model_data.get("submodel", ""),
            "statefreq": statefreq,
            "rates": model_data.get("rates", ""),
            "insertmodel": model_data.get("insertmodel", ""),
            "insertrate": model_data.get("insertrate", "")
        })
        model_counter += 1
    json_data['models'] = models
    
    if 'tree' in user_data:
        json_data['TREE']["template"] = user_data['tree']
    
    if 'b1' in user_data:
        json_data['BRANCHES']['b1'] = user_data['b1']
    
    if 'SEQ_LENGTH' in user_data:
        json_data['PARTITIONS']['part1'][2] = user_data['SEQ_LENGTH']
    
    if 'out_file' in user_data:
        json_data['EVOLVE'][1] = user_data['out_file'].replace("{i}", str(iteration))

    return json_data

# Generate the configuration file
def generate_config(json_data, output_path):
    ###################################################################################################################################
    ######## This function takes the combined json/user_file data and places it into a formatted template for INDELIBLE, ##############
    ######## the output_path is whichever you choose. MAKE SURE YOU KNOW WHERE YOU'RE PLACING IT ######################################
    ###################################################################################################################################

    output_lines = []

    # Add TYPE section
    output_lines.append(f"[TYPE] {json_data['TYPE'][0]}\n\n")

    # Add MODEL sections
    for model in json_data['models']:
        output_lines.append(f"[MODEL] {model['name']}\n")
        if model['submodel']:
            output_lines.append(f"  [submodel] {model['submodel']}\n")
        if model['statefreq']:
            output_lines.append(f"  [statefreq] {model['statefreq']}\n")
        if model['rates']:
            output_lines.append(f"  [rates] {model['rates']}\n")
        if model['insertmodel']:
            output_lines.append(f"  [insertmodel] {model['insertmodel']}\n")
        if model['insertrate']:
            output_lines.append(f"  [insertrate] {model['insertrate']}\n")
        output_lines.append("\n")

    # Add SETTINGS section
    output_lines.append(f"[SETTINGS]\n")
    if 'ancestralprint' in json_data:
        output_lines.append(f"  [ancestralprint] {json_data['ancestralprint'][0]}\n")
    if 'output' in json_data:
        output_lines.append(f"  [output] {json_data['output'][0]}\n")
    if 'randomseed' in json_data:
        output_lines.append(f"  [randomseed] {json_data['randomseed']}\n")
    output_lines.append("\n")

    # Add TREE section
    output_lines.append(f"[TREE] {json_data['TREE']['template']}\n")

    # Add BRANCHES section
    output_lines.append(f"[BRANCHES] b1 {json_data['BRANCHES']['b1']}\n")

    # Add PARTITIONS section
    partitions = json_data['PARTITIONS']['part1']
    output_lines.append(f"[PARTITIONS] part1 [t1 b1 {partitions[2]}]\n")

    # Add EVOLVE section
    output_lines.append(f"[EVOLVE]\npart1 1 {json_data['EVOLVE'][1]}\n")

    # Write to output file
    with open(output_path, 'w') as f:
        for line in output_lines:
            f.write(line)
    print(f"Configuration file {output_path} generated successfully.")
    
###################################################################################################################################
# Paths to the user_file, json template, the dictionaries they go into, and output files                            ###############
##### This is where you should delegate the proper paths, make sure to keep the file names the same, only change ##################
##### the preceeding path names         ###########################################################################################       
###################################################################################################################################


user_file_path = '/home/s36jshem_hpc/sealion/sealion_files/sea_lion/user_file_v2.txt'
json_template_path = '/home/s36jshem_hpc/sealion/sealion_files/INDEL_TEMPLATE_v2.json'

#Read user data from the user_file
user_data, iq_model = read_user_data(user_file_path)


#Load JSON template
with open(json_template_path, 'r') as f:
    json_data = json.load(f)

# Validate user data
validate_user_data(user_data, json_data) ##json_data is the combined user_file into the json template, user_data is the dict with all the user_file info in it

#I added a time stamp to each folder generated, so you always know when you started the run
now = datetime.now()
now_format = now.strftime("%Y-%m-%d_%H-%M-%S")
#I also added the sequence length from the user file to give more details on the script run
BP = user_data['SEQ_LENGTH']
print(BP)

temp_path = f'/home/s36jshem_hpc/sealion/runs/indelible_input_{now_format}_{BP}'
try:
    # exist_ok=True suppresses the exception if folder already exists
    indelible_input_path = temp_path
    os.makedirs(indelible_input_path, exist_ok=True)
except OSError as error:
    print(f"Error! Unable to create directory: {indelible_input_path}/")


# Generate 60 simulation files with different out_file names, 0-10 = 50% GC, 11-20 = 55%... 51-60 = 75% GC
for i in range(1, 61):
    user_data['randomseed'] = random.randint(1, 100000)
    updated_json_data = update_json_structure(json_data, user_data, i)
    output_file_path = f'{indelible_input_path}/indelible_input_{i}.txt'
    generate_config(updated_json_data, output_file_path)

def rename_indel(directory_path, indel_path):
    #################################################################################################################  
    ####This function iterates over all the files you generated earlier, renames them so they can run       #########
    ####in indelible, then renames them back following the run ensuring you don't loose track          ##############
    #################################################################################################################
    
    control_files = []
    for f in os.listdir(directory_path):
        if f.startswith('indelible') and f.endswith('.txt'):
            control_files.append(f)

    for i in control_files:
        original_path = os.path.join(directory_path, i)
        temp_path = os.path.join(indel_path, 'control.txt')
        
        os.rename(original_path, temp_path)
        print('First file processed, renamed to control.txt ')
        
        try:
            #subprocess.run(['./indelible', {temp_path}], cwd = '/home/s36jshem_hpc/sealion/runs' , check=True)
            os.chdir(indel_path)
            command = f"indelible control.txt"
            os.system(command)
        except Exception:
            print("Error loading file properly")
            
        os.rename(temp_path, original_path)
        
        print(f'{i} processed through Indelible, then renamed to original')

indel_path = '/home/s36jshem_hpc/sealion/runs'
directory_path = f'/home/s36jshem_hpc/sealion/runs/indelible_input_{now_format}_{BP}'
rename_indel(directory_path, indel_path)

def move_to_indel_output(model_path, output_directory):
    ##########################################################################################################################
    #### This function moves the output files generated by Indelible to the specified output directory    ################
    ##########################################################################################################################
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    for f in os.listdir(model_path):
        if f.startswith('fasta'):
            src_path = os.path.join(model_path, f)
            dst_path = os.path.join(output_directory, f)
            shutil.move(src_path, dst_path)


indelible_output_path = f'/home/s36jshem_hpc/sealion/runs/indelible_output_{now_format}_{BP}'
model_path  = f'/home/s36jshem_hpc/sealion/runs/'
move_to_indel_output(model_path, indelible_output_path)

def rename_seq_fasta(indelible_output_path, IQTREE_path, iq_model):
     ##########################################################################################################################
     #### This function moves the correct files to iq_output folder, iterates over each file in the iq_output folder, #########
     #### renames the sequence identifiers so they're unique (this is only if your file has repeat sequences ##################
     ##### as they must be unique for IQ-TREE), and then runs each file through IQ-TREE. ######################################
     ##########################################################################################################################  

    indelible_output_path = f'/home/s36jshem_hpc/sealion/runs/indelible_output_{now_format}_{BP}'
    IQTREE_path = f'/home/s36jshem_hpc/sealion/runs/iq_output_{now_format}_{BP}'
    
    # Move files from indelible output to IQ-TREE input
    os.makedirs(IQTREE_path, exist_ok=True)
    pattern = re.compile(r'fastaout\d+\.fas$')
    for f in os.listdir(indelible_output_path):
        if pattern.match(f):
            src_path = os.path.join(indelible_output_path, f)
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
        if f.startswith('fasta') and f.endswith('.fas'):
            indel_out_files.append(f)
    
    for file_name in indel_out_files:
        try:
            file_path = os.path.join(IQTREE_path, file_name)
            command = f"iqtree2 -s {file_path} -m {iq_model} "
            subprocess.run(command, cwd=IQTREE_path, shell=True)
        except Exception as e:
            print(f"Error loading file properly: {e}")


iqtree_output_path = f'/home/s36jshem_hpc/sealion/runs/iq_output_{now_format}_{BP}'
rename_seq_fasta(
    indelible_output_path,
    iqtree_output_path,
    iq_model)

#################################################################################################################################
##The IQTREE outputs 420 files, you have 60 indelible inputs, you then place one fastaout file into the iq_output folder ########
##and it generates 7 files, so 60 * 7 = 420              ########################################################################                                                               #########
##The files are named fastaout1.fas, fastaout2.fas, etc.   ######################################################################                                                             #########
#################################################################################################################################

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

        
move_tree_files(iqtree_output_path, f'/home/s36jshem_hpc/sealion/runs/tree_output_{now_format}_{BP}')

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

    chunk_size = 10
    num_chunks = len(sequences[group])/chunk_size
    
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
    
    shutil.move(clade_path, sealion_directory)
        
    return "Sequences extracted and saved! "

fasta_path = f'/home/s36jshem_hpc/sealion/runs/iq_output_{now_format}_{BP}'
clade_path = f'/home/s36jshem_hpc/sealion/runs/clade_files_{now_format}_{BP}'
sealion_directory = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/'
make_clade_files(fasta_path, clade_path, sealion_directory)

def shrink_newick(newick_string_location):
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
                line1 = re.sub(r':-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', '', line)
                newick_strings_files[file] = line1.strip()
    print(newick_strings_files)

    #### This below replaces the fasta file with the corresponding clade file, so we can run the subprocess command below    
    clade_def = [f for f in os.listdir(f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_{now_format}_{BP}') if f.startswith('clade_def')]

    updated_newick_strings = {}
    for clade_file in clade_def:
        clade_num_match = re.search(r'clade_def_file_(\d+)', clade_file)
        if clade_num_match:
            clade_num = clade_num_match.group(1)
            for fasta_file, newick in newick_strings_files.items():
                fasta_num_match = re.search(r'fastaout(\d+)', fasta_file)
                if fasta_num_match:
                    fasta_num = fasta_num_match.group(1)
                    # Match clade number with fasta number
                    if fasta_num == clade_num:
                        updated_newick_strings[clade_file] = newick
                        #print(f"Matched Clade {clade_file} with Fasta {fasta_file}")
    clade_file_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_{now_format}_{BP}'
    results = {}
    reroot_directory = '/home/s36jshem_hpc/sealion'
    for k, v in updated_newick_strings.items(): ##### All of this below only works SOMETIMES when this script is run in the terminal command line, it works well as an HPC cluster job
        full_clade = os.path.join(clade_file_location, k)
        shutil.move(full_clade, reroot_directory)
        command = f'python3 ESofT.py "{v}" {k}'
        esoft_run = subprocess.run(command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        output = esoft_run.stdout.strip()
        reroot_command = f'./reroot.o "{output}" D {k}'
        reroot_run =  subprocess.run(reroot_command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        reroot_output = reroot_run.stdout.strip()
        results[k] = reroot_output
        dst_path = os.path.join(reroot_directory, k)
        shutil.move(dst_path, clade_file_location)

    sorted_results = dict(sorted(results.items(), key=lambda x: int(re.search(r'clade_def_file_(\d+)', x[0]).group(1))))
    print(sorted_results)
    
    newick_path = f'/home/s36jshem_hpc/sealion/runs/corrected_newick_output_{now_format}_{BP}'
    if not os.path.exists(newick_path):
        os.makedirs(newick_path)

    for clade_def, rerooted_newick in sorted_results.items():
        number_of_file = re.search(r'clade_def_file_(\d+)', clade_def)
        number_for_file = number_of_file.group(1)
        output_file = os.path.join(newick_path, f'corrected_newick_{number_for_file}.txt')
        with open(output_file, 'w') as f:
            f.write(rerooted_newick.strip())
    

tree_location = f'/home/s36jshem_hpc/sealion/runs/tree_output_{now_format}_{BP}'
shrink_newick(tree_location)


def graph_correct_outputs(newick_corrected_path, user_data, tq_dist_path, graph_path):
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
    user_newick = user_data
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = '/home/s36jshem_hpc/sealion/runs/'
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

        command = f"/home/s36jshem_hpc/local/bin/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=tq_dist_path, shell=True, capture_output = True, text = True)
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
    csv_path = os.path.join(graph_path, f'{now_format}.csv')
    savetxt(csv_path, numpy_csv, delimiter=',') 

    #formatted_newick_tree = '\n'.join(textwrap.wrap(user_data['tree'], width=40))
    #formatted_newick_model = '\n'.join(textwrap.wrap(user_data['b1'], width=40))

    #legend_label = f'Newick Branch Lengths:\n{formatted_newick_tree}\n\nNewick Model:\n{formatted_newick_model}'
    #custom_legend = Line2D([0], [0], linestyle = 'None', color='green', label = legend_label)

    plt.plot(x, y, linestyle = '--', label = '_nolegend_')

    #plt.plot(x, regression_line, 'r--', label = 'regression line')

    plt.xlabel('GC Content')
    plt.ylabel('% of Correct Newick Tree Topologies')
    plt.title('Correct Newick String Matches by GC Content (IQ-TREE)')
    plt.xticks(x, gc_content_labels)
    #plt.legend(handles = [custom_legend], loc = 'upper right', fontsize = 'x-small')

    
    plt.savefig(os.path.join(graph_path, f'{now_format}.png'), format='png')
    plt.show()

user_data = "(((A1#F81,B1#F81_2)#F81,C1#F81_2)#F81,D1#F81)#F81;"
newick_path = f'/home/s36jshem_hpc/sealion/runs/corrected_newick_output_{now_format}_{BP}'
tq_dist = '/home/s36jshem_hpc/sealion/runs'
graph_location = f"/home/s36jshem_hpc/sealion/plots/clade_files_{now_format}/"
graph_correct_outputs(newick_path, user_data, tq_dist, graph_location)













#####################################################################
########### Graveyard of coding pasts ###############################     
#####################################################################      
#####      .-"      "-.  ############################################
#####     /            \ ############################################
#####    |              |############################################
#####    |,  .-.  .-.  ,| ###########################################
#####    | )(_o/  \o_)( | ###########################################
#####    |/     /\     \| ###########################################
#####    (_     ^^     _) ###########################################
#####     \__|IIIIII|__/  ###########################################
#####      | \IIIIII/ |   ###########################################
#####      \          /   ###########################################
#####       `--------`  #############################################
##################################################################### 
#####################################################################
#####################################################################

'''
This was supposed to retrieve the log likelihoods
    sub_optimal_newicks = []
    for file in os.listdir(iqtree_output_path):
        if file.endswith('treels'): 
            sub_optimal_newicks.append(file)

    for newick_file in sub_optimal_newicks: 
        likeli_file_match = re.search(r'fastaout(\d+)\.fas\.treels', newick_file)
        if likeli_file_match:
            likeli_file_matched = likeli_file_match.group(1)
            for file_path in indel_out_files:
                og_file_match = re.search(r'fastaout(\d+)\.fas', file_path)
                if og_file_match:
                    og_file_matched = og_file_match.group(1)
                    if og_file_matched == likeli_file_matched:
                        try:
                            command1 = (
                                f"iqtree2 -s {file_path} -m {iq_model} "
                                f"-z {newick_file} -n 0 -redo"
                            )
                            print(f"Now we are running the log-likelihoods of {file_path} with {command1}")
                            subprocess.run(command1, cwd=IQTREE_path, shell=True)
                        except Exception as e:
                            print(f"Error loading log_likelihood file properly: {e}")

'''


'''
def shrink_newick(newick_string_location):
    ###############################################################################################################################################
    #### This function should shrink the newick string down without the branch lengths, then create a dictionary with the fasta file location #####
    ###############################################################################################################################################
    newick_files = [f for f in os.listdir(newick_string_location)]

    full_path = [os.path.join(newick_string_location, f) for f in newick_files]
    newick_strings_files = {}
    for file in full_path:
        if file.endswith('treefile'):
            with open(file, 'r') as f:
                for line in f:
                    line = re.sub(r':[0-9.]+', '', line)
                    newick_strings_files[file] = line.strip()
    
    return newick_strings_files

newick_string_location = f'/home/s36jshem_hpc/sealion/runs/tree_output_{now_format}'
newick_strings_files = shrink_newick(newick_string_location)

def clade_file_sub(clade_file_location):
    #################################################################################################################################
    #### This below replaces the fasta file with the corresponding clade file, so we can run the subprocess command below ###########  
    ################################################################################################################################# 

    clade_def = [f for f in os.listdir(clade_file_location) if f.startswith('clade_def')]
   

    updated_newick_strings = {}
    for clade_file in clade_def:
        clade_num_match = re.search(r'clade_def_file_(\d+)', clade_file)
        if clade_num_match:
            clade_num = clade_num_match.group(1)
            for fasta_file, newick in newick_strings_files.items():
                fasta_num_match = re.search(r'fastaout(\d+)', fasta_file)
                if fasta_num_match:
                    fasta_num = fasta_num_match.group(1)
                    # Match clade number with fasta number
                    if fasta_num == clade_num:
                        updated_newick_strings[clade_file] = newick
                        print(f"Matched Clade {clade_file} with Fasta {fasta_file}")

    return updated_newick_strings

clade_file_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_{now_format}'
updated_newick_strings = clade_file_sub(clade_file_location)

def run_esoft_step(updated_newick_strings, clade_file_location, reroot_directory):
    #################################################################################################################################
    ######This function should run esoft then collect the outputs in a dictionary with the clade file, and the new newick string#####
    #################################################################################################################################
    esoft_outputs = {}
    for clade_file, newick in updated_newick_strings.items():
        full_clade_path = os.path.join(clade_file_location, clade_file)
        try:
            shutil.move(full_clade_path, reroot_directory)
            print(f"Moved {full_clade_path} -> {reroot_directory}")
        except OSError as e:
            if e.errno == errno.EEXIST:
                print(f"File already exists in destination: {clade_file}, skipping.")
                continue
            
        command = f'python3 ESofT.py "{newick}" {clade_file}'
        esoft_run = subprocess.run(command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        print(f'Running: {clade_file} with stripped IQTREE newick: {newick}')
        esoft_output = esoft_run.stdout.strip()
        esoft_outputs[clade_file] = esoft_output
        

        os.chdir(reroot_directory)
        shutil.move(os.path.join(reroot_directory, clade_file), clade_file_location)

    return esoft_outputs

reroot_directory = '/home/s36jshem_hpc/sealion'
esoft_outputs = run_esoft_step(updated_newick_strings, clade_file_location, reroot_directory )

def run_reroot_step(esoft_outputs, reroot_directory, clade_file_location):
    #################################################################################################################################
    ###This function should now run the reroot function, using the outputs from the esoft script which are captured in a dictionary##
    #################################################################################################################################
    rerooted_results = {}

    for clade_file, esoft_output in esoft_outputs.items():
        full_clade_path = os.path.join(clade_file_location, clade_file)
        try:
            shutil.move(full_clade_path, reroot_directory)
            print(f"Moved {full_clade_path} -> {reroot_directory}")
        except OSError as e:
            if e.errno == errno.EEXIST:
                print(f"File already exists in destination: {clade_file}, skipping.")
                continue
            
        reroot_command = f'./reroot.o "{esoft_output}" D {clade_file}'
        
        reroot_run =  subprocess.run(reroot_command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        print(f'Running: {clade_file} with esoft newick: {esoft_output}')

        reroot_output = reroot_run.stdout.strip()
        rerooted_results[clade_file] = reroot_output

        os.chdir(reroot_directory)
        # move clade file back again in case reroot modified it
        clade_path = os.path.join(reroot_directory, clade_file)
        if os.path.exists(clade_path):
            shutil.move(clade_path, clade_file_location)
        if not os.path.exists(clade_path):
            print('Clade path to clade file output incorrectly defined')

    return rerooted_results

rerooted_results = run_reroot_step(esoft_outputs, reroot_directory, clade_file_location)

def write_corrected_newick(rerooted_results, newick_output_directory):
    
    os.makedirs(newick_output_dir, exist_ok=True)

    sorted_results = dict(sorted(rerooted_results.items(), key=lambda x: int(re.search(r'clade_def_file_(\d+)', x[0]).group(1))))
    for clade_def, newick in sorted_results.items():
        match = re.search(r'clade_def_file_(\d+)', clade_def)
        if match:
            file_number = match.group(1)
            output_file = os.path.join(newick_output_directory, f'corrected_newick_{file_number}.txt')
            with open(output_file, 'w') as f:
                f.write(newick.strip())

newick_output_dir = f'/home/s36jshem_hpc/sealion/runs/corrected_newick_output_{now_format}'
write_corrected_newick(rerooted_results, newick_output_dir)
'''

'''
def mv_clade_files(clade_dir, sealion_location): 
    This function moves the files from the clade_path to the sealion location, then runs the sealion script, this is an optionable
    function if you'd like to run each file one at a time and not in parallel.
      
    clade_def_files = [f for f in os.listdir(clade_dir) if f.startswith('clade_def_file') and f.endswith('.txt')]
    clade_files = [f for f in os.listdir(clade_dir) if f.startswith('clade_file') and f.endswith('.fas')]


    clade_def_files.sort(key=lambda x: tuple(map(int, re.search(r'clade_def_file_(\d+)_(\d+)\.txt', x).groups())))
    clade_files.sort(key=lambda x: tuple(map(int, re.search(r'clade_file_(\d+)_(\d+)\.fas', x).groups())))

    if len(clade_def_files) != len(clade_files):
        print("Mismatch between .txt and .fas files")
        return

    for i, (fas_fn, txt_fn) in enumerate(zip(clade_files, clade_def_files), start=1):
        src_fas = os.path.join(clade_dir, fas_fn)
        src_txt = os.path.join(clade_dir, txt_fn)

        try:
            shutil.move(src_fas, sealion_location)
            shutil.move(src_txt, sealion_location)

            os.chdir(sealion_location)
            cmd = (
                f"apptainer exec SeaLion_container_dir sealion1.pl -i {fas_fn} -p {txt_fn} -o D -M '10000' -l '10000' -prt 1 -s"
            )
            print(f"→ Running: {cmd}")
            os.system(cmd)

        except Exception as e:
            print(f"Error processing pair #{i} ({fas_fn}, {txt_fn}): {e}")

    print(f"Processed {i} file pairs through SeaLion")


mv_clade_files(
    f"/home/s36jshem_hpc/sealion/runs/clade_files_{now_format}",
    "/home/s36jshem_hpc/sealion/sealion_script"
)

def make_clade_files(fasta_path, clade_path, sealion_directory):
    this function makes the clade files in a random partitioned manner, randomizing the clade_def file and the clade file for each
    gc step, e.g 1-100, 101-200... then moves them to the sealion portion of the script

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
            output_file = os.path.join(clade_path, f"clade_file_{gc_contents[chunk_idx]}_{j+1}.fas")
            with open(output_file, 'w') as f:
                for group in groups:
                    for label, seq in group_chunks[group][j]:
                        f.write(f">{label}\n{seq}\n")

            clade_def_file = os.path.join(clade_path, f"clade_def_file_{gc_contents[chunk_idx]}_{j+1}.txt")
            with open(clade_def_file, 'w') as f:
                for group in groups:
                    f.write(f"{group}, {', '.join(label_chunks[group][j])}\n")
        
    if not os.path.exists(sealion_directory):
        os.makedirs(sealion_directory)

    shutil.move(clade_path, sealion_directory )

    return "Sequences extracted, saved, and moved!"

fasta_path = f'/home/s36jshem_hpc/sealion/runs/iq_output_{now_format}'
clade_path = f'/home/s36jshem_hpc/sealion/runs/clade_files_{now_format}'
sealion_directory = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir'
make_clade_files(fasta_path, clade_path, sealion_directory)
'''

'''
def run_e_soft(clade_file_location, reroot_directory):
    This function will run the clade file dictionary (clade file and the respective newick string), it will move the clade file into the right directory run the newick
    string, then move the clade file back into the directory and hopefully run that new respective newick string
    clade_file_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_{now_format}'
    results = {}
    reroot_directory = '/home/s36jshem_hpc/sealion'
    for k, v in updated_newick_strings.items():
        full_clade = os.path.join(clade_file_location, k)
        print(full_clade)
        print(reroot_directory)
        shutil.move(full_clade, reroot_directory)
        print(f"Moved {clade_path} -> {reroot_directory}")
        print(k, v)
        command = f'python3 ESofT.py "{v}" {k}'
        esoft_run = subprocess.run(command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        output = esoft_run.stdout.strip()
        print(output)
        reroot_command = f'./reroot.o "{output}" D {k}'
        reroot_run =  subprocess.run(reroot_command, cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        reroot_output = reroot_run.stdout.strip()
        results[k] = reroot_output
        dst_path = os.path.join(reroot_directory, k)
        print(dst_path)
        shutil.move(dst_path, clade_file_location)
    
    sorted_results = dict(sorted(results.items(), key=lambda x: int(re.search(r'clade_def_file_(\d+)', x[0]).group(1))))
    print(sorted_results)


   
    newick_path = f'/home/s36jshem_hpc/sealion/runs/corrected_newick_output_{now_format}'
    if not os.path.exists(newick_path):
        os.makedirs(newick_path)

    for clade_def, rerooted_newick in sorted_results.items():
        number_of_file = re.search(r'clade_def_file_(\d+)', clade_def)
        number_for_file = number_of_file.group(1)
        output_file = os.path.join(newick_path, f'corrected_newick_{number_for_file}.txt')
        with open(output_file, 'w') as f:
            f.write(rerooted_newick.strip())
    

tree_location = f'/home/s36jshem_hpc/sealion/runs/tree_output_{now_format}'
shrink_newick(tree_location)
'''
