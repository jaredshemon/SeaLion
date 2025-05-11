#!/usr/env/bin python3 
import os
import matplotlib
import matplotlib.pyplot as plt
import re
import subprocess
import numpy as np
from datetime import datetime
from numpy import savetxt

now = datetime.now()
now_format = now.strftime("%Y-%m-%d_%H-%M-%S")
 
def diff_visualizations():
    #############################################################################################################################################
    ###This graphs support values on the y axis and the clade file on the x axis, more specifically the most supported topology. The topology ###
    #### is indicated by the color in the legend.                                    ############################################################
    #############################################################################################################################################
    topology_supports = {}
    newick_strings = []

    for j in range(1, 61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-08_13-37-23/sealion_runs/clade_file_{j}.fas/testresults_clade_file_{j}/TSV'
        match = re.search(r'clade_files_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', tsv_location)
        if match:
            clade_file_location = match.group()
        
        best_newick = ''
        best_sup = 0
        
        for data_files in os.listdir(tsv_location):
            if data_files.startswith('MQ1'):
                data_file = data_files
                avg_support_file = os.path.join(tsv_location, data_file)
                with open(avg_support_file, 'r') as f:
                    lines = f.readlines()
                    for i in lines:
                        i.strip()
                        if 'median' in i:
                            list_lines = i.split()
                            sup, newick = list_lines[5], list_lines[1]
                            if float(sup) > best_sup:
                                best_sup = float(sup)
                                best_newick = newick
        if best_newick:
            topology_supports[f'clade_file_{j}'] = [best_sup, best_newick]
            newick_strings.append(best_newick)

    # Build color mapping for unique topologies
    unique_topologies = sorted(set(value[1] for value in topology_supports.values()))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # expand if needed
    color_map = {topo: colors[i % len(colors)] for i, topo in enumerate(unique_topologies)}

    # Prepare plotting
    labels = []
    supports = []
    bar_colors = []

    for clade_file in sorted(topology_supports.keys(), key=lambda x: int(x.split('_')[-1])):
        support, topology = topology_supports[clade_file]
        labels.append(clade_file.replace('clade_file_', ''))
        supports.append(support)
        bar_colors.append(color_map[topology])

    # Plot
    plt.figure(figsize=(14, 6))
    plt.bar(labels, supports, color=bar_colors)

    # Add legend
    legend_handles = [plt.Line2D([0], [0], color=color_map[topo], lw=4, label=topo) for topo in color_map]
    plt.legend(handles=legend_handles, title='Topology')

    plt.xlabel('Clade File #')
    plt.ylabel('Support Value')
    plt.title('Best Supported Topology per Clade File (Colored by Topology)')
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

    saving_location = f"/home/s36jshem_hpc/sealion/plots/{clade_file_location}/"
    if not os.path.exists(saving_location):
        os.makedirs(saving_location)

    plt.savefig(f"{saving_location}/topology_support_plot.png", dpi=300)

    plt.close()

    return best_newick, best_sup, saving_location, newick_strings, clade_file_location

#best_newick, best_sup, saving_location, newick_strings, clade_file_location = diff_visualizations()

def graph_correct_outputs(correct_newick, tq_dist_path):
    ########################################################################################################################################################
    ###This function should take the newick string from our sea_lion output file, cross check it with the original, then graphs the correct ones ###########
    #### vs the incorrect ones.                                    #########################################################################################
    ########################################################################################################################################################
    
    def stripped_newick(string):
        return re.sub(r'([0-9.e-]+|#[A-Za-z0-9_]+)', '', string)

    newick_path = f"{saving_location}"
    newick_file_path = os.path.join(newick_path, 'newickfile1.txt')
    user_newick_path = os.path.join(newick_path, 'newickfile_user.txt' )
    results = []
    

    for i in newick_strings:
        with open(newick_file_path, 'w') as f:
            print(stripped_newick(i))
            f.write(stripped_newick(i))
        with open(user_newick_path, 'w') as f:
            print(stripped_newick(correct_newick))
            f.write(stripped_newick(correct_newick))

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

    #This saves the x,y data as a csv for future graph overlays
    csv_output = {x_axis:y_axis for x_axis,y_axis in zip(list(x),y)}
    numpy_csv = np.array(list(csv_output.items()))
    csv_path = os.path.join(f'{saving_location}', f'Correct_Incorrect.csv')
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

    
    plt.savefig(os.path.join(f'{saving_location}', f'Correct_Incorrect.png'), format='png')
    plt.show()

    return csv_path, clade_file_location

correct_newick = "(((A1#F81,B1#F81_2)#F81,C1#F81_2)#F81,D1#F81)#F81;"
tq_dist = '/home/s36jshem_hpc/sealion/runs'
#csv_path, clade_file_location = graph_correct_outputs(correct_newick, tq_dist)

def overlay_correct(IQ_csv_location):
    ########################################################################################################################################################
    ### This function should overlay the correct vs incorrect from SeaLion v the correct v incorrect from IQTREE ###########################################
    ########################################################################################################################################################
    
    data1 = np.loadtxt(IQ_csv_location, delimiter=',')
    data2 = np.loadtxt(csv_path, delimiter=',')

    x1, y1 = data1[:, 0], data1[:, 1]
    x2, y2 = data2[:, 0], data2[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(x1, y1, marker='o', linestyle='--', color='blue', label='ML(IQTREE)')
    plt.plot(x2, y2, marker='x', linestyle='--', color='green', label='SeaLion')

    plt.xlabel('GC Content Group')
    plt.ylabel('% Correct Topologies')
    plt.title('Overlay of Correct Topology Matches by GC Content')
    plt.xticks(ticks=range(6), labels=[f"{47 + i * 4}%" for i in range(6)])
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(f'{saving_location}/correct_overlay_comparison.png', dpi=300)
    plt.show()

IQ_csv_location = f'/home/s36jshem_hpc/sealion/runs/tree_graphs/2025-05-08_13-37-23/2025-05-08_13-37-23.csv'
#overlay_correct(IQ_csv_location)

def diff_graphs():
    risk_dist_lines = []
    for j in range(1,61):
        tsv_location = f'/home/s36jshem_hpc/sealion/sealion_script/runs_dir/clade_files_2025-05-08_13-37-23/sealion_runs/clade_file_{j}.fas/testresults_clade_file_{j}/TSV'
        for data_files in os.listdir(tsv_location):
            if data_files.startswith('T9') and data_files.endswith('median_nap.tsv'):
                data_file = data_files
                median_file = os.path.join(tsv_location, data_file)
                with open(median_file, 'r') as f:
                    lines = f.readlines()
                    for i in lines:
                        i.strip()
                        if 'RISK_DIST' in i:
                            risk_dist_lines.append(i)
                            print(tsv_location)
                if risk_dist_lines:
                    for r_d_line in risk_dist_lines:
                        print(r_d_line)
diff_graphs()