#!usr/bin/env python3 
# -*- coding: utf-8 -*-

import os
import re
import subprocess
import shutil
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import textwrap
from matplotlib.lines import Line2D
from datetime import datetime
from numpy import savetxt

now = datetime.now()
now_format = now.strftime("%Y-%m-%d_%H-%M-%S")

'''
def shrink_newick(newick_string_location):
    This function should shrink the newick string down without the branch lengths, then it will run it through juliannes python script 'ESOFT' then it 
    will run it through reroot
    newick_files = [f for f in os.listdir(newick_string_location)]
    
    full_path = [os.path.join(tree_location, f) for f in newick_files]
    
    newick_strings_files = {}
    for file in full_path:
        with open(file, 'r') as f:
            for line in f:
                line = re.sub(r':[0-9.]+', '', line)
                newick_strings_files[file] = line.strip()

    #### This below replaces the fasta file with the corresponding clade file, so we can run the subprocess command below    
    clade_def = [f for f in os.listdir('/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-05_19-48-29') if f.startswith('clade_def')]

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
    clade_file_location = '/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-05_19-48-29'
    results = {}
    reroot_directory = '/home/s36jshem_hpc/sealion'
    for k, v in updated_newick_strings.items():
        full_clade = os.path.join(clade_file_location, k)
        shutil.move(full_clade, reroot_directory)
        command = f'python3 ESofT.py "{v}" {k}'
        esoft_run = subprocess.run([command], cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        output = esoft_run.stdout.strip()
        reroot_command = f'./reroot.o "{output}" D {k}'
        reroot_run =  subprocess.run([reroot_command], cwd = reroot_directory, capture_output = True, check=True, text = True, shell = True)
        reroot_output = reroot_run.stdout.strip()
        results[k] = reroot_output
        dst_path = os.path.join(reroot_directory, k)
        shutil.move(k, clade_file_location)

    sorted_results = dict(sorted(results.items(), key=lambda x: int(re.search(r'clade_def_file_(\d+)', x[0]).group(1))))
    print(sorted_results)
    
    newick_path = '/home/s36jshem_hpc/sealion/runs/corrected_newick_output_2025-05-05_19-48-29'
    if not os.path.exists(newick_path):
        os.makedirs(newick_path)

    for clade_def, rerooted_newick in sorted_results.items():
        number_of_file = re.search(r'clade_def_file_(\d+)', clade_def)
        number_for_file = number_of_file.group(1)
        output_file = os.path.join(newick_path, f'corrected_newick_{number_for_file}.txt')
        with open(output_file, 'w') as f:
            f.write(rerooted_newick.strip())
    

tree_location = '/home/s36jshem_hpc/sealion/runs/tree_output_2025-05-05_19-48-29'
shrink_newick(tree_location)
'''

def graph_correct_outputs(newick_corrected_path, user_data, tq_dist_path, graph_path):
    ''' This function should take the newick string from our tree file, cross check it with the original, then graphs the correct ones
    vs the incorrect ones.'''
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
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
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
    incorrect_counts = [results[i:i + 10].count(1) for i in range(0, 60, 10)]
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
    plt.title('Correct Newick String Matches by GC Content')
    plt.xticks(x, gc_content_labels)
    #plt.legend(handles = [custom_legend], loc = 'upper right', fontsize = 'x-small')

    
    plt.savefig(os.path.join(graph_path, f'{now_format}.png'), format='png')
    plt.show()

user_data = "(((A1#F81,B1#F81_2)#F81,C1#F81_2)#F81,D1#F81)#F81;"
newick_path = f'/home/s36jshem_hpc/sealion/runs/corrected_newick_output_2025-05-06_10-12-51'
tq_dist = '/home/s36jshem_hpc/sealion/runs'
graph_location = f"/home/s36jshem_hpc/sealion/runs/tree_graphs/{now_format}"
graph_correct_outputs(newick_path, user_data, tq_dist, graph_location)