# Phylogenetic Simulation and Topology Evaluation Pipeline

## Overview

This repository contains a comprehensive, end-to-end computational pipeline developed for a Master's thesis in Organismic Biology. The pipeline is designed to **simulate sequence alignments under controlled GC-content heterogeneity**, infer phylogenetic trees using **IQ-TREE**, evaluate competing topologies using **SeaLion**, and quantitatively compare inferred trees against a known ground-truth topology using **quartet distance metrics**. The workflow is optimized for execution on an **HPC environment** and emphasizes reproducibility, scalability, and structured output suitable for downstream statistical and graphical analyses.

The pipeline integrates simulation, inference, topology filtering, and visualization into a single orchestrated script while maintaining modularity through well-defined functional stages.

---

## Scientific Purpose

This pipeline was developed to address the following research questions:

* How does **GC-content variation** across clades affect phylogenetic inference accuracy?
* How do **filtered vs. unfiltered** topology support metrics differ under heterogeneous evolutionary scenarios?
* How does **IQ-TREE likelihood support** compare to **SeaLion RISK + DIST filtering** when identifying the correct topology?

Ground-truth trees are known *a priori*, allowing explicit classification of inferred topologies as **correct** or **incorrect**.

---

## Pipeline Summary

The workflow proceeds through the following major stages:

1. **Alignment Simulation (AliSim via IQ-TREE)**
2. **Tree Inference (IQ-TREE)**
3. **Clade File Construction**
4. **Topology Filtering and Support Estimation (SeaLion)**
5. **Tree Re-rooting and Canonicalization (ESOFT + reroot)**
6. **Quartet Distance Evaluation**
7. **Statistical Summaries and Visualization**

Each stage writes structured intermediate outputs to disk, enabling checkpointing and independent re-analysis.

---

## Dependencies

### Required Software

The following external tools must be installed and accessible in the execution environment:

* **Python ≥ 3.8**
* **IQ-TREE 2** (compiled with AliSim support)
* **SeaLion** (executed via Apptainer/Singularity container)
* **ESOFT** (tree canonicalization script)
* **reroot** (binary for outgroup re-rooting)
* **quartet_dist / tq_dist** (quartet distance calculation)

### Python Libraries

The pipeline relies on standard scientific Python packages:

* `numpy`
* `pandas`
* `matplotlib`
* `multiprocessing`
* `statistics`

All Python dependencies are compatible with Conda-based environments commonly used on HPC clusters.

---

## Input Configuration

### Command-Line Arguments

The script expects the **working directory** as its first positional argument:

```bash
python3 pipeline.py <working_directory>
```

### User-Defined Paths

Key external dependencies are specified at the top of the script, including:

* SeaLion container location
* IQ-TREE substitution model
* Location of reroot binaries
* Path to quartet distance executable
* AliSim user configuration file

These must be edited prior to execution.

### AliSim Configuration File

Simulation parameters (sequence length, GC ranges, substitution model, invariant sites, etc.) are read from a structured text file. The pipeline dynamically expands GC-content ranges and embeds them into Newick templates for simulation.

---

## Output Structure

Each execution generates timestamped directories to ensure reproducibility:

```
working_directory/
├── ALI_output_<timestamp>/          # Simulated alignments and Newick trees
├── iq_output_<timestamp>/           # IQ-TREE inference outputs
├── tree_output_<timestamp>/         # Inferred treefiles
├── corrected_IQ_newick_output_/     # Canonicalized and rerooted trees
├── runs_dir/
│   └── clade_files_<timestamp>/     # SeaLion input/output directories
├── plots/                           # All figures and summary CSVs
└── log_folder_<timestamp>/          # Execution logs
```

All plots are exported in **SVG format** and all numerical summaries are saved as **CSV files** for downstream analysis.

---

## Key Analyses and Visualizations

The pipeline produces multiple classes of figures and tables, including:

* **Percentage of Correct Trees vs. GC Content (IQ-TREE)**
* **Best Supported Topology per Dataset (Filtered and Unfiltered SeaLion)**
* **Support Threshold Categorization (Good / Moderate / High Conflict)**
* **Δ Support Between Best and Second-Best Topologies**
* **Comparison of Correct vs. Incorrect Topologies**
* **IQ-TREE Log-Likelihood Δ Analyses (log-scaled)**

Each figure corresponds to a numbered table and SVG output suitable for thesis figures.

---

## Parallelization

SeaLion analyses are executed using Python’s `multiprocessing` module, with one process per dataset (default: 60). This design is optimized for HPC cluster execution and significantly reduces wall-clock time.

---

## Error Handling and Logging

* All runtime errors are captured via a custom exception handler.
* Detailed logs are written to a timestamped log directory.
* Failed subprocess calls are reported with sufficient context for debugging.

This design ensures traceability of failures during long-running cluster jobs.

---

## Reproducibility Notes

* All randomization is confined to AliSim simulations.
* Ground-truth Newick strings are explicitly defined and preserved.
* Outputs are deterministic given identical inputs and software versions.

For archival purposes, it is recommended to store large simulation outputs externally (e.g., institutional storage or Zenodo) rather than in Git repositories.

---

## Intended Audience

This repository is intended for:

* Phylogeneticists studying model misspecification and compositional bias
* Computational biologists conducting large-scale simulation studies
* Researchers evaluating topology filtering and support metrics

Familiarity with phylogenetic inference tools and HPC environments is assumed.

---

## Author

**Jared Shemonsky**
Master’s Student, Organismic Biology
University of Bonn

---
