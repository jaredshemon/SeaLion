#!/usr/env/bin python3 
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




now = datetime.now()
now_format = now.strftime("%Y-%m-%d_%H-%M-%S")
 
#############################################################################################################################################
###This graphs support values on the y axis and the clade file on the x axis, more specifically the most supported topology. The topology ###
#### is indicated by the color in the legend.                                    ############################################################
############################################################################################################################################# 
def diff_visualizations():
    
    topology_supports = {}
    unfiltered_topology_supports = {}
    newick_strings = [] #FILTERED NEWICKS FROM SEALION (if fully rejected == 0)
    newick_strings1 = [] #UNFILTERED NEWICKS FROM SEALION

    for j in range(1, 61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-15_20-59-58_10000/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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

    # Plot
    plt.figure(figsize=(14, 6))
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
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    # Add legend
    legend_handles = [plt.Line2D([0], [0], color=color_map[topo], lw=4, label=topo) for topo in color_map if topo != "N/A"]
    asterisk_handle = mlines.Line2D([], [], color='none', marker='*', markersize=14, 
                                markerfacecolor='red', label='Unsupported', linestyle='None')
    legend_handles.append(asterisk_handle)
    plt.legend(handles=legend_handles, bbox_to_anchor=(1, 0.95))

    plt.xlabel('Dataset')
    plt.ylabel('Support Value (Max 1)')
    plt.title('Best Supported RISK+DIST Filtered SeaLion Topology per Clade File (Colored by Topology)')
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

    saving_location = f"/home/s36jshem_hpc/sealion/plots/{clade_file_location}/"
    if not os.path.exists(saving_location):
        os.makedirs(saving_location)

    plt.savefig(f"{saving_location}/Best_Supported_RISK+DIST_Topology", dpi=300)

    plt.close()

    return best_newick, best_sup, saving_location, newick_strings1, clade_file_location, clade_file_time, tsv_location, unfiltered_topology_supports, newick_strings

best_newick, best_sup, saving_location, newick_strings1, clade_file_location, clade_file_time, tsv_location, unfiltered_topology_supports, newick_strings = diff_visualizations()

################################################################################
###This should graph the same as above albeit the UNFILTERED supports  #########
################################################################################  
def unfiltered_quartet_supports(unfiltered_topology_supports):

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

    # Plot
    plt.figure(figsize=(14, 6))
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
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    # Add legend
    legend_handles = [plt.Line2D([0], [0], color=color_map[topo], lw=4, label=topo) for topo in color_map]
    plt.legend(handles=legend_handles, bbox_to_anchor=(1, 0.95))

    plt.xlabel('Dataset')
    plt.ylabel('Support Value (Max 1)')
    plt.title('Best Supported Unfiltered SeaLion Topology per Clade File (Colored by Topology)')
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

    saving_location = f"/home/s36jshem_hpc/sealion/plots/{clade_file_location}/"
    if not os.path.exists(saving_location):
        os.makedirs(saving_location)

    plt.savefig(f"{saving_location}/Best_Supported_Unfiltered_Topology", dpi=300)

    plt.close()

unfiltered_quartet_supports(unfiltered_topology_supports)

###############################################################################################
###This should graph the same topology graph as the two above but for the IQTREE analyis  #####
###############################################################################################
def IQ_quartet_supports():

    IQ_likeli_loc = f"/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-15_20-59-58_10000"
    for file in os.listdir(IQ_likeli_loc):
        if file.startswith('fastaout') and file.endswith('ckp.gz'):
            os.chdir(IQ_likeli_loc)
            os.system(f'gunzip {file}')

    newick_scores = []
    diff_newick_scores = {}
    for file in os.listdir(IQ_likeli_loc):
        if file.startswith('fastaout') and file.endswith('ckp'): 
            file = os.path.join(IQ_likeli_loc, file)
            with open(file, 'r') as f:
                lines = f.read()
                log_likelihoods = re.findall(r'\d+:\s+(-\d+\.\d+)', lines)
                top = sorted(log_likelihoods)[0]
                diff_newick_scores[file] = [top]

    y_labels = []
    x_labels = []
    for file_name, score in diff_newick_scores.items():
        match = re.search(r'fastaout(\d+)\.fas\.ckp', file_name)
        if match:
            dataset_number = int(match.group(1))
            x_labels.append(dataset_number)
            y_labels.append(float(score[0]))

    
    newick_corrected_path = '/home/s36jshem_hpc/sealion/runs/corrected_newick_output_2025-05-15_20-59-58_10000'
    user_newick = '(((A,B),C),D);'
    tq_dist_path = '/home/s36jshem_hpc/local/bin/quartet_dist'
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

    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = '/home/s36jshem_hpc/sealion/runs/'
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    
    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            f.write(i)
        with open(user_newick_path, 'w') as f:
            f.write(user_newick)

        command = f"/home/s36jshem_hpc/local/bin/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=newick_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))
    
    sorted_data = sorted(zip(y_labels, x_labels))
    x_labels_sorted, y_values_sorted = zip(*sorted_data)

    x = np.arange(1,61)
    # Assign colors based on results (0 = correct, 1 = incorrect)
    colors = ['1f77b4' if result == 0 else 'blue' for result in results]

    plt.figure(figsize=(16, 6))
    plt.bar(x, y_labels, color=colors, align='center')

    correct_newick = set(newick_strings[i] for i, result in enumerate(results) if result == 0)
    incorrect_newick = set(newick_strings[i] for i, result in enumerate(results) if result == 1)

    correct_legend_label = "".join(correct_newick)
    incorrect_legend_label = "".join(incorrect_newick)

    correct_patch = plt.Line2D([0], [0], color='orange', lw=4, label=correct_legend_label)
    incorrect_patch = plt.Line2D([0], [0], color='blue', lw=4, label=incorrect_legend_label)
    plt.legend(handles=[correct_patch, incorrect_patch], loc='upper right')

    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, -77900, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)


    plt.xlabel("Dataset Number")
    plt.ylabel("Log-Likelihood Score")
    plt.ylim(-75000,-78000)
    plt.title("Log-Likelihood Scores by Dataset")
    plt.xticks(range(1,61))
    plt.tight_layout()
    plt.show()

    plt.savefig(f"{saving_location}/IQ_Topology_barchart", dpi=300)

    plt.close()

#IQ_quartet_supports()

###############################################################################################################################################################
###This function should take the UNFILTERED newick string from our SeaLion output file, cross check it with the original, then graphs the correct ones ########
#### vs the incorrect ones.                                    ################################################################################################
###############################################################################################################################################################
def graph_correct_outputs(newick_strings1, correct_newick, tq_dist_path):

    
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = f"{saving_location}"
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    

    for i in newick_strings1: #Newick strings contains the unfiltered score newick strings, and comes from the diff_visualizations() function
        with open(newick_file_path, 'w') as f:
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            f.write(stripped_newick(correct_newick))


        command = f"/home/s36jshem_hpc/local/bin/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=tq_dist_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))

    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
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
    csv_path = os.path.join(f'{saving_location}/Tree_Success_GC_Content_SeaLion.csv')
    savetxt(csv_path, numpy_csv, delimiter=',') 

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

    
    plt.savefig((f'{saving_location}/Tree_Success_GC_Content_SeaLion.png'))
    plt.show()


    return csv_path, clade_file_location, results

correct_newick = "(((A,B),C),D);"
tq_dist = '/home/s36jshem_hpc/sealion/runs'
csv_path, clade_file_location, results = graph_correct_outputs(newick_strings1, correct_newick, tq_dist)

###############################################################################################################################################################
###This function should take the FILTERED newick string from our SeaLion output file, cross check it with the original, then graphs the correct ones ##########
#### vs the incorrect ones.                                    ################################################################################################
###############################################################################################################################################################
def graph_correct_outputs2(newick_strings, correct_newick, tq_dist_path):

    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = f"{saving_location}"
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

        
            command = f"/home/s36jshem_hpc/local/bin/quartet_dist {newick_file_path} {user_newick_path}"
            result = subprocess.run(command, cwd=tq_dist_path, shell=True, capture_output = True, text = True)
            output = result.stdout.strip()
            results_filtered.append(int(output))
            rejected_focused.append('unrejected')
        else:
            results_filtered.append(int(1))
            rejected_focused.append('rejected')

    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    correct_counts = [results_filtered[i:i + 10].count(0) for i in range(0, 60, 10)]
    incorrect_counts = [results_filtered[i:i + 10].count(1) for i in range(0, 60, 10)]
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
    csv_path = os.path.join(f'{saving_location}/RISKDIST_Tree_Success_GC_Content_SeaLion.csv')
    savetxt(csv_path, numpy_csv, delimiter=',') 

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

    
    plt.savefig((f'{saving_location}/RISKDIST_Tree_Success_GC_Content_SeaLion.png'))
    plt.show()


    return csv_path, clade_file_location, results_filtered, rejected_focused

correct_newick = "(((A,B),C),D);"
tq_dist = '/home/s36jshem_hpc/sealion/runs'
#csv_path1, clade_file_location, results_filtered, rejected_focused = graph_correct_outputs2(newick_strings, correct_newick, tq_dist)

#############################################################################################################################################
###This function should take the newick string from our tree file, cross check it with the original, then graphs the correct ones ###########
#### vs the incorrect ones. (IQTREE)                              ###########################################################################
#############################################################################################################################################    
def graph_correct_outputsIQ(newick_corrected_path, correct_newick_string_user_data, tq_dist_path):
    
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
    
 
    # Preparing the data for graphing
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
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
    
    plt.savefig(f'{saving_location}/IQTREE_Success_GC_Content.png', format='png')
    plt.show()

corrected_newick_path = f'/home/s36jshem_hpc/sealion/runs/corrected_newick_output_2025-05-15_20-59-58_10000'
#graph_correct_outputsIQ(corrected_newick_path, correct_newick, tq_dist)

#################################################################
### This function just replots the IQTREE correct quartets  #####                                              
#################################################################
def IQ_correct(IQ_csv_location):

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
    plt.savefig(f'{saving_location}/Tree_Success_GC_Content_IQTREE.png', dpi=300)
    plt.show()

    return x1, y1

#IQ_csv_location = f'/home/s36jshem_hpc/sealion/plots/clade_files_2025-05-15_20-59-58/2025-05-15_20-59-58.csv'
#IQ_correct(IQ_csv_location)

##################################################################################################################
#### This function should overlay the correct vs incorrect vs rejected from Unfiltered SeaLion          ##########
##################################################################################################################
def correct_incorrect_rejected_filtered(results_filtered, rejected_focused):
    gc_contents = [47 + (i // 10) * 4 for i in range(60)]
    gc_content_labels = [f"{47 + i * 4}%" for i in range(6)]
    correct_counts = [results_filtered[i:i + 10].count(0) for i in range(0, 60, 10)]
    correct_percent = [int(i)/10 for i in correct_counts]
    incorrect_counts = [results_filtered[i:i + 10].count(1) for i in range(0, 60, 10)]
    incorrect_percent = [int(i)/10 for i in incorrect_counts]
    percent_counts = [(correct / (correct + incorrect)) if (correct + incorrect) > 0 else 0 for correct, incorrect in zip(correct_counts, incorrect_counts)]
    rejected_counts = [0 if i == 'unrejected' else 1 for i in rejected_focused]
    rejected_counts1 = [rejected_counts[i:i+10].count(1) for i in range(0,60,10)]
    rejected_percent= [int(i)/10 for i in rejected_counts1]

    x = range(6)
    plt.figure(figsize=(10, 6))
    plt.plot(x, correct_percent, linestyle='--', marker='o', color='green', label='Correct Topologies')
    plt.plot(x, incorrect_percent, linestyle='--', marker='o', color='darkgoldenrod', label='Incorrect Topologies')
    plt.plot(x, rejected_percent, linestyle='--', marker='o', color='darkviolet', label='Rejected Topologies')

    plt.xlabel('GC Content')
    plt.ylabel('% Tree Success')
    plt.ylim(0,1)
    plt.title('Tree Success, Rejection, and Failure (SeaLion)')
    plt.xticks(x, gc_content_labels)
    plt.legend()

    
    plt.savefig((f'{saving_location}correct_incorrect_rejected_SeaLion.png'))
    plt.show()

#correct_incorrect_rejected_filtered(results_filtered, rejected_focused)

###################################################################################################################################
#### This function should overlay the correct vs incorrect from Unfiltered SeaLion v the correct v incorrect from IQTREE ##########
###################################################################################################################################
def overlay_correct(csv_path, IQ_csv_location):
    
    data1 = np.loadtxt(IQ_csv_location, delimiter=',')
    data2 = np.loadtxt(csv_path, delimiter=',')


    x1, y1 = data1[:, 0], data1[:, 1]
    x2, y2 = data2[:, 0], data2[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='ML(IQTREE)')
    plt.plot(x2, y2, marker='o', linestyle='--', color='green', label='SeaLion')

    plt.xlabel('GC Content Group')
    plt.ylabel('% Correct Topologies')
    plt.title('Overlay of Correct Topology Matches by GC Content (Unfiltered)')
    plt.xticks(ticks=range(6), labels=[f"{47 + i * 4}%" for i in range(6)])
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Unfiltered_Overlay_Tree_Success.png', dpi=300)
    plt.show()

    return x1, y1, x2, y2

#IQ_csv_location = f'/home/s36jshem_hpc/sealion/plots/clade_files_2025-05-15_20-59-58/IQTREE_Success_GC_Content.csv'
#x1, y1, x2, y2 = overlay_correct(csv_path, IQ_csv_location)

###################################################################################################################################
#### This function should overlay the correct vs incorrect from Filtered SeaLion v the correct v incorrect from IQTREE ############
###################################################################################################################################
def overlay_correct2(csv_path1, IQ_csv_location):
    
    data1 = np.loadtxt(IQ_csv_location, delimiter=',')
    data2 = np.loadtxt(csv_path1, delimiter=',')


    x1, y1 = data1[:, 0], data1[:, 1]
    x2, y2 = data2[:, 0], data2[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='ML(IQTREE)')
    plt.plot(x2, y2, marker='o', linestyle='--', color='green', label='SeaLion')

    plt.xlabel('GC Content Group')
    plt.ylabel('% Correct Topologies')
    plt.title('Overlay of Correct Topology Matches by GC Content (RISK+DIST)')
    plt.xticks(ticks=range(6), labels=[f"{47 + i * 4}%" for i in range(6)])
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Overlay_Tree_Success.png', dpi=300)
    plt.show()

    return x1, y1, x2, y2

#x1, y1, x2, y2 = overlay_correct2(csv_path1, IQ_csv_location)

#################################################################################################################################################
### This graphs the difference between the best and second best tree topologys from SeaLion ONLY POSITIVE BECAUSE THE DICT IS SORTED ############
#################################################################################################################################################
def diff_graphs(tsv_location):

    diff = []
    diff1 = []
    for j in range(1,61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-15_20-59-58_10000/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
    
    plt.figure(figsize=(16, 5))
    plt.plot(y_axis, differences, marker='+', linestyle='--', color='red', label='Sealion Support Δ' )
    incomplete_x = [i+1 for i,v in enumerate(differences) if v == 0]
    incomplete_y = [0] * len(incomplete_x)
    plt.scatter(incomplete_x, incomplete_y, marker='*', color='black', s=180, label='Fully Rejected')
    plt.axhline(y=.2, color='blue', linestyle=':', linewidth=1, label=f'High Conflict Below')


        # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    plt.xlabel('Dataset')
    plt.ylabel('Support Δ')
    plt.title('Sealion Δ of Best Topology v. Second Best Topology (RISK+DIST)')
    #plt.xticks(ticks=range(50), labels=[f"{47 + i * 4}%" for i in range(5)])
    plt.xticks(ticks=range(1, 61))
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/SeaLion_best_second_Δ.png', dpi=300)
    plt.show()

    return differences, differencesU
    
differences, differencesU = diff_graphs(tsv_location)

############################################################################################
### Same as the graph above but for unfiltered data ONLY POSITIVE ##########################
############################################################################################
def diff_graphs1(differencesU):

    y_axis = range(1, 61)
    plt.figure(figsize=(16, 5))
    plt.plot(y_axis, differencesU, marker='+', linestyle='--', color='red', label='Sealion Support Δ' )
    plt.axhline(y=.2, color='blue', linestyle=':', linewidth=1, label=f'High Conflict Below')
        # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='lightcoral', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .99, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    plt.xlabel('Dataset')
    plt.ylabel('Support Δ')
    plt.title('Sealion Δ of Best Topology v. Second Best Topology (Unfiltered)')
    #plt.xticks(ticks=range(50), labels=[f"{47 + i * 4}%" for i in range(5)])
    plt.xticks(ticks=range(1, 61))
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/SeaLion_Unfil_best_second_Δ.png', dpi=300)
    plt.show()

#diff_graphs1(differencesU)

############################################################################################
### Shows the delta when the tree is correct v. incorrect ONLY POSITIVE ####################
############################################################################################
def diff_tree_correct_v_incorrect(differencesU, results):
    for j in range(1,61):
        newicks_scores = {}
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-15_20-59-58_10000/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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

    plt.figure(figsize=(18, 5))
    x_axis = range(1,61)
    y_axis = [float(i) for i in differencesU]
    for x, res  in zip(x_axis, results):
        if res == 0:
            plt.bar(x, y_axis[x-1], color = 'skyblue', label = 'Correct Topology Δ' )
        else:
            plt.bar(x, y_axis[x-1], color = 'seagreen', label = 'Incorrect Topology Δ')


    # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .95, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
    threshold = .2
    plt.axhline(y=threshold, color='lightcoral', linestyle=':', linewidth=1)
    plt.text(61, 0.1, 'High Conflict', color='black', fontsize=10, va='top')


    plt.xlabel('Dataset Index')
    plt.ylabel('Δ of Support Values')
    plt.title('Δ Correct Topology v. Incorrect Topology')
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.ylim(0,1)
    legend_patch = mpatches.Patch(color='skyblue', label='Correct Topology Δ')
    legend_patch1 = (mpatches.Patch(color='seagreen', label='Incorrect Topology Δ'))
    plt.legend(handles=[legend_patch, legend_patch1])  
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Δ_Correct_Topology_v._Incorrect_Topology.png', dpi=300)
    plt.show()

diff_tree_correct_v_incorrect(differencesU, results)
########################################################################################################################################################
### This graphs the difference between the best and second best tree topologys from IQTREE indicated by log-likelihood        ##########################
########################################################################################################################################################
def diff_graphs2(IQ_likeli_loc):

    IQ_likeli_loc = f"/home/s36jshem_hpc/sealion/runs/iq_output_2025-05-15_20-59-58_10000"
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
    plt.plot(indices, diffs, marker='+', linestyle='--', color='red', label='IQTREE ΔLnl ' )

    # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='red', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    #GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, 1e-4, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)
    
    # Highlight significant differences
    threshold = 1e-1
    significant = diffs >= threshold
    plt.axhline(y=threshold, color='lightcoral', linestyle=':', linewidth=1, label=f'Δ ≥ {threshold}')

    plt.xlabel('Dataset')
    plt.ylabel('Δ Log-Likelihood (Best - 2nd Best)')
    plt.title('ΔLnl (IQ-TREE)')
    plt.xticks(indices, rotation=90)
    plt.yscale('log')  # log scale to emphasize small differences
    plt.ylim(bottom=1e-5, top=max(diffs)*1.5)
    plt.legend()
    #plt.tight_layout()
    plt.savefig(f'{saving_location}/IQ_best_second_Δ.png', dpi=300)
    plt.show()

    return diffs, indices

#IQ_likeli_loc = f'/home/s36jshem_hpc/sealion/runs/iq_output_{clade_file_time}'
#diffs, indices = diff_graphs2(IQ_likeli_loc)

##########################################################################################
### Combined graph with two y-axes to compare filtered and unfiltered datasets ###########
##########################################################################################
def combined_graph(differences, differencesU):

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
    plt.axhline(y=.2, color='mediumvioletred', linestyle=':', linewidth=1, label=f'High Conflict Below')

    # Add dashed lines and shaded regions for filtered data
    for x in range(10, 60, 10):
        ax1.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)
    for i in range(0, 60, 20):  # Shaded regions
        ax1.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    # GC content annotations
    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        ax1.text(i * 10 + 5, 0.98, f'{gc}%', ha='center', va='top', fontsize=9, transform=ax1.transData)

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
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Combined_SeaLion_Delta_Support.png', dpi=300)
    plt.show()

# Call the function with your data
#combined_graph(differences, differencesU)

#################################################################################################
### Combined graph with two y-axes to compare filtered and unfiltered datasets Barchart #########           
#################################################################################################
def combined_graph_bar(differences, differencesU):

    x_axis = np.arange(1, 61)  # Shared x-axis for both datasets
    bar_width = 0.4  # Width of each bar

    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(16, 6))

    # First dataset (filtered data) on the left y-axis
    ax1.bar(x_axis - bar_width / 2, differences, width=bar_width, color='darkgoldenrod', label='Filtered Δ Support')
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Support Δ (Filtered)', color='darkgoldenrod')
    ax1.tick_params(axis='y', labelcolor='darkgoldenrod')
    plt.axhline(y=.2, color='green', linestyle=':', linewidth=1, label=f'Good Support')
    plt.text(61, 0.1, 'High Conflict', color='black', fontsize=10, va='top')

    for x in range(1, 60, 10):
        plt.axvline(x=x - 0.5, color='darkgoldenrod', linestyle='--', linewidth=1)

    for i in range(1, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    gc_labels = [47, 51, 55, 59, 63, 67]
    for i, gc in enumerate(gc_labels):
        plt.text(i * 10 + 5, .98, f'{gc}%', ha='center', va='top', fontsize=9, transform=plt.gca().transData)

    # Plot the unfiltered dataset (differencesU) as bars on the right y-axis
    ax2 = ax1.twinx()
    ax2.bar(x_axis + bar_width / 2, differencesU, width=bar_width, color='seagreen', label='Unfiltered Δ Support')
    ax2.set_ylabel('Support Δ (Unfiltered)', color='seagreen')
    ax2.tick_params(axis='y', labelcolor='seagreen')

    # Add legends for both datasets
    fig.legend(loc='upper right', bbox_to_anchor=(.98, .9), bbox_transform=ax1.transAxes)

    # Title and layout
    plt.title('Comparison of Filtered and Unfiltered Δ Support')
    plt.xticks(ticks=range(1, 61))
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Barchart_Combined_SeaLion_Delta_Support.png', dpi=300)
    plt.show()

# Call the function with your data
#combined_graph_bar(differences, differencesU)

############################################################################################
### Combined graph with two y-axes to compare filtered and IQ_likelihood delta datasets  ###
############################################################################################
def combined_graph_IQ(differences, diffs, indices, saving_location):

    y_axis = range(1, 61)  # Shared x-axis for both datasets
    gc_labels = [47, 51, 55, 59, 63, 67]

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
    plt.savefig(f'{saving_location}/Combined_Delta_SeaLion_LogLikelihood.png', dpi=300)
    plt.show()

#combined_graph_IQ(differences, diffs, indices, saving_location)

##############################################################################################
### Combined graph with two y-axes to compare unfiltered and IQ_likelihood delta datasets ####
##############################################################################################
def combined_graph_IQ_Unfil(differencesU, diffs, indices, saving_location):

    y_axis = range(1, 61)  # Shared x-axis for both datasets
    gc_labels = [47, 51, 55, 59, 63, 67]

    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(16, 6))

    # Plot filtered dataset on the left y-axis
    ax1.plot(y_axis, differencesU, linestyle='--', color='seagreen', label='Filtered Δ Support')
    ax1.set_xlabel('Dataset')
    ax1.set_ylabel('Support Δ (Filtered)', color='seagreen')
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
    plt.savefig(f'{saving_location}/Combined_Support_Delta_Unfiltered_LogLikelihood.png', dpi=300)
    plt.show()

#combined_graph_IQ_Unfil(differencesU, diffs, indices, saving_location)

################################################################################################################
### This should look at the difference between correct topology support before and after filtering Barchart ####        
################################################################################################################
def support_b4_af_filtering(results):
    after_support = []
    before_support = []
    rejected = []
    for j in range(1,61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-15_20-59-58_10000/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
    colors = ['darkorange' if i > 0 else 'skyblue' for i in y_1]
    plt.figure(figsize=(12, 6))

    # Plot before and after support
    plt.bar(x_axis, y_1, color=colors, edgecolor='black', label='Difference before v. after filtering')
    label_added = False  # Flag to add the label only once


    for x, rej in zip(x_axis, rejected):
        if rej == 1:
            plt.scatter(x, 0, marker="x", color="black", s=100, label="rejected" if not label_added else "")
            label_added = True
    
    for x, res in zip(x_axis, results):
        if res == 1:
            plt.bar(x, y_1[x-1], color='red', edgecolor='black', label='Failed Topologies')
    
    x_handle = mlines.Line2D([], [], color='none', marker='x', markersize=10, 
                        markerfacecolor='red', label='Rejected', linestyle='None')

    # Labels and title
    plt.xlabel('Dataset Index')
    plt.ylabel('Support Value')
    plt.title('Δ Correct Topology After Filtering vs. Before Filtering')
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.ylim(-0.05, 0.12)
    legend_patch = mpatches.Patch(color='gray', label='Difference After v. Before filtering')
    legend_patch1 = (mpatches.Patch(color='red', label='Failed Topologies'))
    plt.legend(handles=[legend_patch, x_handle, legend_patch1])  
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/Correct_Topology_b4_After_Filtering_Bar.png', dpi=300)
    plt.show()

support_b4_af_filtering(results)

##########################################################################################################
### This should look at the correct topology support before and after filtering LineGraph (DONT LOVE) ####     
##########################################################################################################
def support_b4_af():

    after_support = []
    before_support = []
    for j in range(1,61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-15_20-59-58_10000/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
    plt.savefig(f'{saving_location}/SeaLion_correct_topology_b4_after_filtering.png', dpi=300)
    plt.show()
                           
support_b4_af()

#### THIS DOESNT WORK RIGHT!!! TRY AND EDIT IT FOR EACH INDIVIDUAL DATASET 
#######################################################################
### This graphs the rejected trees as a function of the GC contents ###   
#######################################################################
def reject_GC():

    clade_file_rejected = defaultdict(list)
    percent_rejected = []
    accepted = []
    for j in range(1,61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-15_20-59-58_10000/sealion_runs/clade_file_{j}.fas/testresult_clade_file_{j}/TSV'
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
                                percent_rejected.append(rejected1)
                                clade_file_rejected[clade_file_location] = float(rejected1)
                            elif 'initially':
                                accepted1 = parts[4]
                                accepted.append(accepted1)
   
    print(percent_rejected)
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
    plt.bar(x_axis, y_axis, color='seagreen', edgecolor='black', label='Perecent of Topologies Rejected')


    # Add dashed lines
    for x in range(10, 60, 10):
        plt.axvline(x=x - 0.5, color='blue', linestyle='--', linewidth=1)

    #Add shaded backgrounds
    for i in range(0, 60, 20):  # every other bin
        plt.axvspan(i - 0.5, i + 9.5, color='gray', alpha=0.1)

    # Add labels, title, and legend
    plt.xlabel('Dataset')
    plt.ylabel('% Rejected')
    plt.title('% of SeaLions Rejected Topologies vs GC Content')
    plt.ylim(0, 1)
    plt.xticks(ticks=list(x_axis), labels=[str(i) for i in x_axis], rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Save and show
    plt.savefig(f'{saving_location}/Bar_SeaLion_percent_rejected.png', dpi=300)
    plt.show()

reject_GC()










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
        for scores in identifier_scores.values():
            top_two = sorted(scores, reverse = True)[0:2]
            diff.append(top_two[0] - top_two[1])
        
        if diff:
            avg.append(sum(diff) / len(diff))
'''

'''
    newick_corrected_path = '/home/s36jshem_hpc/sealion/runs/corrected_newick_output_2025-05-15_20-59-58_10000'
    user_newick = '(((A,B),C),D);'
    tq_dist_path = '/home/s36jshem_hpc/local/bin/quartet_dist'
    #def extract_number(file_name):
        #match = re.search(r'\\d+', file_name)
        #return int(match.group()) if match else -1

    newick_strings = []
    for file in sorted(os.listdir(newick_corrected_path), key = extract_number):
        if file.endswith('txt'):
            file_path = os.path.join(newick_corrected_path, file)
            with open(file_path, 'r') as f:
                newick_string = f.read().strip()
                newick_strings.append(newick_string)

    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = '/home/s36jshem_hpc/sealion/runs/'
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    
    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            f.write(i)
        with open(user_newick_path, 'w') as f:
            f.write(user_newick)

        command = f"/home/s36jshem_hpc/local/bin/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=newick_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))
    
    print(results)
'''

'''

    newick_corrected_path = '/home/s36jshem_hpc/sealion/runs/corrected_newick_output_2025-05-15_20-59-58_10000'
    user_newick = '(((A,B),C),D);'
    tq_dist_path = '/home/s36jshem_hpc/local/bin/quartet_dist'
    def extract_number(file_name):
        match = re.search(r'\\d+', file_name)
        return int(match.group()) if match else -1

    newick_strings = []
    for file in sorted(os.listdir(newick_corrected_path), key = extract_number):
        if file.endswith('txt'):
            file_path = os.path.join(newick_corrected_path, file)
            with open(file_path, 'r') as f:
                newick_string = f.read().strip()
                newick_strings.append(newick_string)

    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = '/home/s36jshem_hpc/sealion/runs/'
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    
    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            f.write(i)
        with open(user_newick_path, 'w') as f:
            f.write(user_newick)

        command = f"/home/s36jshem_hpc/local/bin/quartet_dist {newick_file_path} {user_newick_path}"
        result = subprocess.run(command, cwd=newick_path, shell=True, capture_output = True, text = True)
        output = result.stdout.strip()
        results.append(int(output))
    
    y = range(len(results))
    x = results
    plt.figure(figsize=(14, 6))
    plt.bar(y, x, width =.8)
    
    plt.xlabel('Dataset')
    plt.ylabel('Correct v. Incorrect')
    plt.title('Correct v. Incorrect IQTREE Datasets')
    plt.xticks(rotation=90)
    plt.ylim(1, 0)
    plt.tight_layout()
    plt.show()

    plt.savefig(f"{saving_location}/IQ_Correct_v_Incorrect_barchart", dpi=300)

    plt.close()
'''