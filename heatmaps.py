
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


def percent_rejected_heatmaps():
#######################################################################################
### This graphs a heatmap of the percent rejected trees for each dataset in trs1-4 ####   
#######################################################################################

    clade_file_rejected = defaultdict(list)
    percent_rejected = []
    accepted = []
    now_format = ['2025-06-19_10-53-07', '2025-06-17_15-02-16', '2025-06-19_10-56-42', '2025-06-17_15-02-43', '2025-06-17_15-14-13', '2025-06-17_15-53-24']
    for i in range(1,7):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1'
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
    for i in range(6):
        start = i * 60
        end = start + 60
        trs.append([
            list(percent_rejected[j].values())[0]
            for j in range(start, end)
        ])
    trs1, trs2, trs3, trs4, Trs5, Trs6 = trs 
    print(trs1, trs2, trs3, trs4, Trs5, Trs6)
    
    fig = plt.figure(figsize=(30,30))
    ax = sns.heatmap(
        trs,
        cmap='coolwarm',
        annot=False,
        fmt=".1f",
        cbar_kws={'label': 'Percent Rejected'},
        xticklabels=np.arange(1, 61),
        yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6']
    )

    # Make the colorbar label bigger
    cbar = ax.collections[0].colorbar
    cbar.set_label('Percent Rejected', fontsize=35)   
    cbar.ax.tick_params(labelsize=35)                 
    ax.set_title("Percent Rejected for Each Dataset (Trs1–6)", fontsize = 35, fontweight = 'bold')
    ax.set_xlabel("Dataset", fontsize = 35)
    ax.set_ylabel("Topology (Trs1–6)", fontsize = 35)
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6' ], rotation=0, fontsize = 35)

    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='x', which='both', top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(5, 60, 10))   # position ticks at 10, 20, ..., 60
    ax.set_xticklabels(np.arange(1, 7), fontsize=30, fontweight='bold')
    for x in np.arange(10, 51, 10):
        ax.axvline(x, color='black', linestyle='-', linewidth=5, alpha=1)

    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/percent_rejected_heatmap_trs1-6.svg', dpi=400)
        
percent_rejected_heatmaps()

def delta_support_heatmaps():
##############################################################################################################
### This should look at the difference between correct topology support before and after filtering heatmap ###   
##############################################################################################################
    rejected1 = []
    differences = []
    now_format = ['2025-06-19_10-53-07', '2025-06-17_15-02-16', '2025-06-19_10-56-42', '2025-06-17_15-02-43', '2025-06-17_15-14-13', '2025-06-17_15-53-24']
    for i in range(1,7):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1'
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
    for i in range(6):
        start = i * 60
        end = start + 60
        trs.append([
            list(differences[j].values())[0]
            for j in range(start, end)
        ])
    trs1, trs2, trs3, trs4, trs5, trs6 = trs 
    
    rejected_mask = np.array([list(d.values())[0] for d in rejected1])  # shape (4, 60)
    masked_trs = np.where(rejected_mask == 1, np.nan, trs)

    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(
        trs,
        cmap='coolwarm',
        annot=False,
        fmt=".1f",
        cbar_kws={'label': 'Support Δ'},
        xticklabels=np.arange(1, 61),
        yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6']
    )

    # Increase colorbar label and tick font size
    cbar = ax.collections[0].colorbar
    cbar.set_label('Support Δ', fontsize=20)   
    cbar.ax.tick_params(labelsize=16)          
    ax.set_title("Δ Correct Topology After Filtering vs. Before Filtering (Trs1–6)", fontsize = 20, fontweight = 'bold')
    ax.set_xlabel("Dataset", fontsize = 15)
    ax.set_ylabel("Topology (Trs1–6)", fontsize = 18)
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6' ], rotation=0, fontsize = 18)

    for i in range(rejected_mask.shape[0]):  # 4 trs
        for j in range(rejected_mask.shape[1]):  # 60 datasets
            if rejected_mask[i, j] == 1:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color='black', linewidth=0))
    rejected_patch = mpatches.Patch(color='black', label='Rejected')
    ax.legend(
        handles=[rejected_patch], loc='upper left', bbox_to_anchor=(1.2, 1), borderaxespad=0, frameon=False, ncol=1)    
    cbar = ax.collections[0].colorbar

            
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/delta_b4_after_support_trs1-6.svg', dpi=400)
    
#delta_support_heatmaps()

def unfiltered_sup_heatmaps():
#######################################################################################################
###This graphs a heatmap of the support values for unfiltered, and collects the values for filtered ###
#######################################################################################################
    now_format = ['2025-06-19_10-53-07', '2025-06-17_15-02-16', '2025-06-19_10-56-42', '2025-06-17_15-02-43', '2025-06-17_15-14-13', '2025-06-17_15-53-24']
    topology_supports= []
    unfiltered_topology_supports = []

    for k in range(1,7):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{k}.1'
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
    for i in range(6):
        start = i * 60
        end = start + 60
        trs.append([
            list(unfiltered_topology_supports[j].values())[0]
            for j in range(start, end)
        ])
    trs1, trs2, trs3, trs4, trs5, trs6 = trs 

    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='coolwarm', annot=False, fmt=".1f", cbar_kws={'label': 'Unfiltered Support Values'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'])
    ax.set_title("Unfiltered Support Values for Each Dataset (Trs1–6)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–6)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'], rotation=0)
    print("minimum", np.min(trs) )
    
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/unfiltered_heatmap_trs1-6.svg', dpi=400)
    
    return topology_supports, unfiltered_topology_supports

#unfiltered_topology_supports, topology_supports = unfiltered_sup_heatmaps()

def filtered_sup_heatmaps(topology_supports):
    print(topology_supports)
    trs = []  
    for i in range(6):
        start = i * 60
        end = start + 60
        trs.append([
            list(topology_supports[j].values())[0]
            for j in range(start, end)
        ])
    
    trs1, trs2, trs3, trs4, trs5, trs6 = trs 
    fig = plt.figure(figsize=(20,20))
    ax = sns.heatmap(trs, cmap='coolwarm', annot=False, fmt=".1f", cbar_kws={'label': 'Filtered Support Values'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'])
    ax.set_title("Filtered Support Values for Each Dataset (Trs1–6)")
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Topology (Trs1–6)")
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'], rotation=0)
    print("minimum", np.min(trs) )
    rejected_mask = np.full((6, 60), np.nan)
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
    plt.savefig('/home/jshemonsky/sealion/heatmaps/filtered_heatmap_trs1-6.svg', dpi=400)
    
    return topology_supports

#filtered_sup_heatmaps(topology_supports)

def IQ_support_heatmap():
    ###############################################################################################
    ###This should graph the same topology graph as the two above but for the IQTREE analyis  #####
    ###############################################################################################
    user_newick = "(((A,B),C),D);"
    tq_dist_path = "/home/jshemonsky/sealion/tqdist/tqDist-1.0.2/install/bin/"
    sorted_data1 = []
    for i in range(1,7):
        now_format = ['2025-06-19_10-53-07', '2025-06-17_15-02-16', '2025-06-19_10-56-42', '2025-06-17_15-02-43', '2025-06-17_15-14-13', '2025-06-17_15-53-24']
        newick_path = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1'
        IQ_likeli_loc = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1/iq_output_{now_format[i-1]}'
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
        newick_corrected_path = f'/home/jshemonsky/sealion/runs/runs/trs{i}.1/corrected_IQ_newick_output_{now_format[i-1]}'
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
    trs1, trs2, trs3, trs4, trs5, trs6 = [
        [list(sorted_data1[k].keys())[j] for j in range(60)]
        for k in range(6)
    ]
    trs = [trs1, trs2, trs3, trs4, trs5, trs6]    
    
    fig = plt.figure(figsize=(30,30))
    ax = sns.heatmap(
        trs,
        cmap='bwr',
        annot=False,
        fmt=".1f",
        cbar_kws={'label': 'IQTREE Support Values'},
        xticklabels=np.arange(1, 61),
        yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6']
    )

    # Increase colorbar label + tick font size
    cbar = ax.collections[0].colorbar
    cbar.set_label('IQTREE Support Values', fontsize=35)   
    cbar.ax.tick_params(labelsize=35)    
    ax.set_title("IQTREE Support Values for Each Dataset (Trs1–6)", fontsize = 35, fontweight = 'bold')
    ax.set_xlabel("Dataset", fontsize = 35)
    ax.set_ylabel("Topology (Trs1–6)", fontsize = 35)
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'], rotation=0, fontsize = 35)
    
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='x', which='both', top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(5, 60, 10))   # position ticks at 10, 20, ..., 60
    ax.set_xticklabels(np.arange(1, 7), fontsize=30, fontweight='bold')
    for x in np.arange(10, 51, 10):
        ax.axvline(x, color='black', linestyle='-', linewidth=5, alpha=1)
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/IQ_supports_heatmap_trs1-6.svg', dpi=400)
    
IQ_support_heatmap()

def IQ_delta_sup():
    ########################################################################################################################################################
    ### This graphs the difference between the best and second best tree topologys from IQTREE indicated by log-likelihood        ##########################
    ########################################################################################################################################################    
    sorted_data1 = []
    diff_newick_scores = []
    for i in range(1,7):
        now_format = ['2025-06-19_10-53-07', '2025-06-17_15-02-16', '2025-06-19_10-56-42', '2025-06-17_15-02-43', '2025-06-17_15-14-13', '2025-06-17_15-53-24']
        newick_path = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1'
        IQ_likeli_loc = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1/iq_output_{now_format[i-1]}'
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
    trs = [diff_newick_scores[i:i + 60] for i in range(0, 360, 60)]
    trs = [[list(d.values())[0] for d in group] for group in trs]
    trs1, trs2, trs3, trs4, trs5, trs6 = trs 

    fig = plt.figure(figsize=(22,22))

    ax = sns.heatmap(
        trs,
        cmap='plasma',
        norm=LogNorm(),
        annot=False,
        fmt=".1f",
        cbar_kws={'label': 'IQTREE Δ Support Values'},
        xticklabels=np.arange(1, 61),
        yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6']
    )

    # Adjust colorbar label and tick font size
    cbar = ax.collections[0].colorbar
    cbar.set_label('IQTREE Δ Support Values', fontsize=30)    
    cbar.ax.tick_params(labelsize=30) 
    ax.set_title("IQTREE Δ Support Values for Each Dataset (Trs1–6)", fontsize = 30, fontweight = 'bold')
    ax.set_xlabel("Dataset", fontsize = 30)
    ax.set_ylabel("Topology (Trs1–6)", fontsize = 30)
    ax.set_yticklabels(['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'], rotation=0, fontsize = 30)
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='x', which='both', top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(5, 60, 10))   # position ticks at 10, 20, ..., 60
    ax.set_xticklabels(np.arange(1, 7), fontsize=30, fontweight='bold')
    for x in np.arange(10, 51, 10):
        ax.axvline(x, color='black', linestyle='-', linewidth=5, alpha=1)
    plt.tight_layout()
    
    plt.tight_layout()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/IQ_delta_supports_heatmap_trs1-6.svg', dpi=400)
    
IQ_delta_sup()

def delta_sup_v_tot_sup_heatmap(unfiltered_topology_supports):
###########################################################################################################
### This the total support of the newick tree v. the delta support of the best v. the second best #########   
###########################################################################################################
        
    now_format = ['2025-06-19_10-53-07', '2025-06-17_15-02-16', '2025-06-19_10-56-42', '2025-06-17_15-02-43', '2025-06-17_15-14-13', '2025-06-17_15-53-24']
    diff = []
    diff1 = []
    for i in range(1,7):
        folder_location = f'/home/jshemonsky/sealion/runs/runs/tbl.01_ibl.05/trs{i}.1'
 
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
    trs5 = []
    trs6 = []

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
        elif u_keys and d_keys == 5:
            trs5.append({u_values:d_values}) 
        elif u_keys and d_keys == 6:
            trs6.append({u_values:d_values}) 
                    
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

    datasets = [trs1, trs2, trs3, trs4, trs5, trs6]
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
    datasets = [trs1, trs2, trs3, trs4, trs5, trs6]  # Each is a list of dicts

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
    ax = sns.heatmap(region_index_grid, cmap=cmap, annot=False, fmt=".1f", cbar_kws={'label': 'Support/Difference Regions by Custom Color'},    xticklabels=np.arange(1, 61), yticklabels=['Trs1', 'Trs2', 'Trs3', 'Trs4', 'Trs5', 'Trs6'])

    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=region_colors[lab], label=lab) for lab in region_labels]
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.15, 1), loc='upper left', borderaxespad=0., ncol = 1)

    plt.tight_layout()
    plt.show()
    plt.savefig('/home/jshemonsky/sealion/heatmaps/delta_support_v_tot_sup_unfiltered_heatmap_trs1-6.svg', dpi=400)
       
#delta_sup_v_tot_sup_heatmap(unfiltered_topology_supports)

