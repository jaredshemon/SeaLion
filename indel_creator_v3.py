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

now_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

#This below is for the logging folder
log_folder = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/log_folder_{now_time}'
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
    '''This reads the userfile information into a dictionary, along with the models which need to be handled
    separately so that they don't overwrite eachother'''
   
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
    '''This funciion validates that whatever is in the ancestralprint, output, and model_name parameters is in 
    the json template'''
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
    '''This function places the dictionary: user_data's information into the json template, the model function at the 
    bottom is needed so the user can input minGC, maxGC, and stepGC, and get a steadily increasing GC content according 
    to the min/max/step every 100 iterations'''
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
                if (iteration - 1) // 100 == (i - gc_min) // gc_step:
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
    '''This function takes the combined json/user_file data and places it into a formatted template for INDELIBLE,
    the output_path is whichever you choose. MAKE SURE YOU KNOW WHERE YOU'RE PLACING IT'''

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
    output_lines.append(f"[TREE] {json_data['TREE']["template"]}\n")

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

# Paths to the user_file, json template, the dictionaries they go into, and output files
'''This is where you should delegate the proper paths, make sure to keep the file names the same, only change
    the preceeding path names'''
    

user_file_path = '/Users/jaredshemonsky/Documents/Thesis Documents/user_file_v2.txt'
json_template_path = '/Users/jaredshemonsky/Documents/Thesis Documents/INDEL_TEMPLATE_v2.json'

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

temp_path = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/indelible_input_{now_format}'
try:
    # exist_ok=True suppresses the exception if folder already exists
    indelible_input_path = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/indelible_input_{now_format}'
    os.makedirs(indelible_input_path, exist_ok=True)
except OSError as error:
    print(f"Error! Unable to create directory: {indelible_input_path}/")


# Generate 600 simulation files with different out_file names, 0-100 = 50% GC, 101-200 = 55%... 501-600 = 75% GC
for i in range(1, 601):
    user_data['randomseed'] = random.randint(1, 100000)
    updated_json_data = update_json_structure(json_data, user_data, i)
    output_file_path = f'{indelible_input_path}/indelible_input_{i}.txt'
    generate_config(updated_json_data, output_file_path)

def rename_indel(directory_path):
    '''This function iterates over all the files you generated earlier, renames them so they can run 
    in indelible, then renames them back following the run ensuring you don't loose track'''
    import os
    import subprocess
    
    control_files = []
    for f in os.listdir(directory_path):
        if f.startswith('indelible') and f.endswith('.txt'):
            control_files.append(f)

    for i in control_files:
        original_path = os.path.join(directory_path, i)
        temp_path = os.path.join(directory_path, 'control.txt')
        os.rename(original_path, temp_path)
        print('First file processed, renamed to control.txt ')
        
        try:
            subprocess.run(["/Users/jaredshemonsky/Downloads/INDELibleV1.03/src/indelible", "control.txt"], cwd = directory_path, shell = True)
        except Exception:
            print("Error loading file properly")
            
        os.rename(temp_path, original_path)
        
        print(f'{i} processed through Indelible, then renamed to original')
        
        
rename_indel(indelible_input_path)


def move_to_indel_output(directory_path, output_directory):
    '''This function moves the output files generated by Indelible to the specified output directory.'''
 
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    for f in os.listdir(directory_path):
        if f.startswith('GTR'):
            src_path = os.path.join(directory_path, f)
            dst_path = os.path.join(output_directory, f)
            shutil.move(src_path, dst_path)


indelible_output_path = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/indelible_output_{now_format}'
move_to_indel_output(indelible_input_path, indelible_output_path)

def rename_seq_fasta(indelible_output_path, IQTREE_path, IQTREE_exe_path, iq_model):
    '''This function moves the correct files to iq_output folder, iterates over each file in the iq_output folder, 
       renames the sequence identifiers so they're unique (this is only if your file has repeat sequences 
       as they must be unique for IQ-TREE), and then runs each file through IQ-TREE.'''
        
    indelible_output_path = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/indelible_output_{now_format}'
    IQTREE_path = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/iq_output_{now_format}'
    IQTREE_exe_path = '/Users/jaredshemonsky/Downloads/iqtree-2.4.0-macOS/bin/iqtree2'

    # Move files from indelible output to IQ-TREE input
    os.makedirs(IQTREE_path, exist_ok=True)
    pattern = re.compile(r'GTRout\d+\.fas$')
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

    # Uncomment and fix the following block if you need to rename sequence identifiers
    '''
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

    # Run IQ-TREE on each file in the IQ-TREE input directory
    indel_out_files = []
    for f in os.listdir(IQTREE_path):
        if f.startswith('GTR') and f.endswith('.fas'):
            indel_out_files.append(f)
    
    for file_name in indel_out_files:
        try:
            file_path = os.path.join(IQTREE_path, file_name)
            command = f"{IQTREE_exe_path} -s {file_path} -m {iq_model}"
            subprocess.run(command, cwd=IQTREE_path, shell=True)
        except Exception as e:
            print(f"Error loading file properly: {e}")


iqtree_output_path = f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/iq_output_{now_format}'
rename_seq_fasta(
    indelible_output_path,
    iqtree_output_path,
    '/Users/jaredshemonsky/Downloads/iqtree-2.4.0-macOS/bin/iqtree2',
    iq_model)

##The IQTREE outputs 4200 files, you have 600 indelible inputs, you then place one GTRout file into the iq_output folder
##and it generates 7 files, so 600 * 7 = 4200
##The files are named GTRout1.fas, GTRout2.fas, etc.

def move_tree_files(iqtree_output_path, tree_output_path):
    '''This function moves the tree files generated by IQ-TREE to the specified tree_output directory.'''
    
    if not os.path.exists(tree_output_path):
        os.makedirs(tree_output_path)   
    
    for f in os.listdir(iqtree_output_path):
        if f.endswith('.treefile'):
            src_path = os.path.join(iqtree_output_path, f)
            dst_path = os.path.join(tree_output_path, f)
            shutil.move(src_path, dst_path)

        
move_tree_files(iqtree_output_path, f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/tree_output_{now_format}')


def graph_correct_outputs(tree_output_path, user_data, tq_dist_path, graph_path):
    ''' This function should take the newick string from our tree file, cross check it with the original, then graphs the correct ones
    vs the incorrect ones.'''
    def extract_number(file_name):
        match = re.search(r'\d+', file_name)
        return int(match.group()) if match else -1


    newick_strings = []
    for file in sorted(os.listdir(tree_output_path), key = extract_number):
        if file.endswith('treefile'):
            file_path = os.path.join(tree_output_path, file)
            with open(file_path, 'r') as f:
                newick_string = f.read().strip()
                newick_strings.append(newick_string)
    user_newick = user_data['b1']
    
    def stripped_newick(string):
        return re.sub(r'(:[0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = '/Users/jaredshemonsky/Downloads/tqDist-1.0.2'
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []


    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            f.write(stripped_newick(user_newick))

        command = f"quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=tq_dist_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))

 
    # Preparing the data for graphing
    gc_contents = [47 + (i // 100) * 4 for i in range(600)]
    gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    correct_counts = [results[i:i + 100].count(0) for i in range(0, 600, 100)]
    incorrect_counts = [results[i:i + 100].count(1) for i in range(0, 600, 100)]
    percent_counts = [(correct / (correct + incorrect)) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]

    
    # Plotting the results
    x = range(6)
    y = percent_counts
    coefficients = np.polyfit(x, y, deg=1)
    polynomial = np.poly1d(coefficients)
    regression_line = polynomial(x)

    formatted_newick_tree = '\n'.join(textwrap.wrap(user_data['tree'], width=40))
    formatted_newick_model = '\n'.join(textwrap.wrap(user_data['b1'], width=40))

    legend_label = f'Newick Branch Lengths:\n{formatted_newick_tree}\n\nNewick Model:\n{formatted_newick_model}'
    custom_legend = Line2D([0], [0], linestyle = 'None', color='green', label = legend_label)

    plt.plot(x, y, linestyle = '--', label = '_nolegend_')

    #plt.plot(x, regression_line, 'r--', label = 'regression line')

    plt.xlabel('GC Content')
    plt.ylabel('% of Correct Newick Tree Topologies')
    plt.title('Correct Newick String Matches by GC Content')
    plt.xticks(x, gc_content_labels)
    plt.legend(handles = [custom_legend], loc = 'upper right', fontsize = 'x-small')

    if not os.path.exists(graph_path):
        os.makedirs(graph_path)
    plt.savefig(os.path.join(graph_path, f'{now_format}.png'), format='png')
    plt.show()


graph_correct_outputs(f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/tree_output_{now_format}', user_data, '/Users/jaredshemonsky/Downloads/tqDist-1.0.2', f'/Users/jaredshemonsky/Downloads/INDELibleV1.03/tree_graphs/{user_data['tree']}')




