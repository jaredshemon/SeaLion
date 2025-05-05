#!usr/bin/env python3 
# -*- coding: utf-8 -*-

import os
import re
import subprocess
import shutil

def shrink_newick(newick_string_location):
    '''This function should shrink the newick string down without the branch lengths, then it will run it through juliannes python script 'ESOFT' then it 
    will run it through reroot'''
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