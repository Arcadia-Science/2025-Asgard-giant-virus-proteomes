**2025-Asgard-giant-virus-proteomes**

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This repository contains the scripts and analysis notebooks used to construct and analyze a deeply annotated  proteome database for Asgard archaea and Giant Viruses (GV, or NCLDV). The primary motivation is to leverage integrated functional, topological, localization, and structural context predictions, alongside homology searches, to:

Characterize the boundaries of sequence divergence with respect to conservation of structure and function.

Identify and prioritize structurally uncharacterized proteins for future study.

Facilitate novel discovery into protein sequence/structure/function relationships

Identify orthologous groups relevant to the evolutionary relationships between Asgard archaea and eukaryotes (detailed phylogenetic analysis is part of a parallel project).

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda here. After installing conda and mamba, run the following command to create the primary pipeline run environment used for most analyses. Note: Specific tools like DeepTMHMM and USPNet required separate environments/setups (Docker, specific Python versions).

# Environment used for data processing, analysis notebooks, MMseqs2 etc.

mamba env create -n asgard_gv_env --file envs/dev.yml
conda activate asgard_gv_env

(Note: The envs/dev.yml file will need to be created based on the packages used, e.g., pandas, numpy, matplotlib, seaborn, biopython, mmseqs2, awscli, etc.)

Install your pre-commit hooks:

pre-commit install

This installs the pre-commit hooks defined in your config (./.pre-commit-config.yaml).

Export your conda environment before sharing:

As your project develops, the number of dependencies in your environment may increase. Whenever you install new dependencies (using either pip install or mamba install), you should update the environment file using the following command.

conda env export --no-builds > envs/dev.yml
```--no-builds` removes build specification from the exported packages to increase portability between different platforms.

## Data

## Input Data:

Source Proteomes: FASTA files for Asgard archaea (311 genomes), Giant Viruses (451 genomes), Eukaryotes (63 genomes), TACK/Euryarchaeota (outgroups). (Specify source/accessions or deposition if applicable).

Reference Data: interpro_entry.list (from InterPro), UniProtKB/RefSeq data used implicitly by tools, AFDB_seq_db, PDB_seq_db

Mapping Files: mapping_parquet_proteinid_to_uniprotkb_or_upi.tsv.

## Intermediate Data:

Filtered FASTA files (e.g., Asgard_all_globular_proteins.fasta).

OrthoFinder output directories (e.g., data/Asgard_Orthofinder_Results, data/GV_Orthofinder_Results).

InterProScan output directories (TSV, GFF3, e.g., InterProScan_Results/Asgard_annotated_nr90/, InterProScan_Results/unknowns_scan/).

USPNet intermediate/output directories (eg. USPNet_Processed_Data*/results.csv).

AlphaFold DB search files: afdb_found_uniprot_acs_or_upi.csv

MMseqs2 database indices (e.g., MGnify_DB/).

MMseqs2 search results (results_vs_mgnify.m8).

PDB sequence homology Search Results: results_vs_pdb_v2.m8.

Extracted sequence/ID lists (e.g., unique_virus_names.txt, afesm_esm_only_uniprot_ids.txt).

## Primary Output Data:

proteome_database_vX.Y.csv: Iterative versions of the main integrated database. Latest version contains comprehensive annotations.

all_filtered_out_proteins_vX.Y.csv: Processed versions of proteins initially filtered based on PDB/AFDB hits.

output_plots/: Directory containing generated figures from analysis notebooks.

## Data Deposition: (TODO: Add Zenodo DOI or other repository links if data is deposited).

Zenodo (tbd)

## Overview

### Description of the folder structure

.
├── .git/                     # Git internal tracking directory (hidden)
├── .gitignore                # Files and directories ignored by Git
├── README.md                 # This file: project description and overview
├── envs/                     # Conda environment definition files
│   └── dev.yml               # Environment specification for reproducibility
├── notebooks/                # Jupyter notebooks for analysis and visualization
│   ├── 00_Setup_Environment.ipynb # Setup (colors, fonts)
│   ├── 01_Data_Exploration.ipynb  # Initial analysis (Composition, Length)
│   ├── 02_Annotation_Analysis.ipynb # Func Cat, IPR, SP/Localization analysis
│   ├── 03_Unknown_Analysis.ipynb  # Investigation of Unknown/Domain categories
│   └── 04_Divergence_Analysis.ipynb # Length variation, Domain Arch comparison
│   └── ...                   # Additional notebooks
├── scripts/                  # Standalone Python or Shell scripts
│   ├── add_specific_category_IPR_v10.py # Custom function categorization
│   ├── extract_fasta.py      # Script to create FASTA files from CSV
│   ├── run_uspnet.sh         # Wrapper script for USPNet
│   ├── run_deeptmhmm.sh      # Wrapper script for DeepTMHMM
│   ├── add_disorder.py       # Script to run Metapredict and add column
│   └── ...                   # Other processing scripts
├── data/                     # Input data (only smaller reference files in Git)
│   └── reference/            # Reference files needed by scripts/notebooks
│       ├── interpro_entry.txt
│       ├── mapping_parquet_proteinid_to_uniprotkb_or_upi.tsv
│       └── afdb_found_uniprot_acs_or_upi.csv
│       └── ...               # Other small (< few MB) lookup tables
│   └── *(Note: Large input FASTAs/DBs stored externally, e.g., S3, EBI FTP)*
└── results/                  # Output files (generally NOT committed to Git - add to .gitignore)
    ├── tables/               # Output data tables (e.g., final DB CSVs)
    │   ├── proteome_database_v1.1.csv
    │   └── all_filtered_out_proteins_v0.8.csv # Example latest version
    │   └── ...
    ├── plots/                # Generated figures and plots
    │   └── output_plots_arcadia_style/ # Directory created by setup notebook
    │       ├── setup_test_plot.png
    │       └── ...
    ├── search_results/       # Homology search output files
    │   ├── results_vs_pdb_v2.m8
    │   └── results_vs_mgnify.m8
    └── predictions/          # Outputs from prediction tools
        ├── DeepTMHMM_Results_Retry/
        └── USPNet_Processed_Data*/

### Methods

Data Acquisition & Preparation: Downloaded source proteomes (Asgard: 311, GV: 451). Filtered proteomes for length (80-1000 aa) and predicted globularity (Metapredict). Generated custom FASTA headers.

Orthology Inference: Ran OrthoFinder (via tack_analysis_env Docker image) separately on the filtered Asgard and GV proteomes to assign proteins to Orthogroups (OGs).

Functional Annotation (InterProScan):

Ran InterProScan (v5.73-104.0 via Docker) iteratively. Initially on nr90 subsets (created via cd-hit -c 0.90), then on proteins initially lacking hits to maximize coverage. Used a curated application list (Pfam, CDD, Gene3D, SUPERFAMILY, SMART, ProSitePatterns, ProSiteProfiles) and requested GO terms/pathways. Managed memory via modified interproscan.sh.

Custom Functional Categorization: Developed and applied a rule-based Python script (add_specific_category_IPR_v10.py) using specific IPR IDs, keywords (in IPR names/source annotations), and IPR type fallback to assign Specific_Functional_Category and Category_Trigger.

Taxonomy Refinement: Fetched NCBI TaxIDs via Entrez; parsed Asgard Phylum, Virus Name; assigned Virus Family via keyword matching.

Sequence Feature Prediction:

Ran USPNet locally (Python env tack_env_x86 / asgard_gv_env) to predict signal peptides (Signal_Peptide_USPNet, SP_Cleavage_Site_USPNet).

Ran DeepTMHMM via Docker on AWS EC2 to predict transmembrane topology (Status: Needs re-run/debug).

Ran Metapredict locally (Python env asgard_gv_env) via add_disorder.py script to calculate Percent_Disorder.

Database Integration: Merged all annotations and metadata into CSV files using Python/pandas scripts. Generated Mature_Protein_Sequence, Predicted_Subcellular_Localization, Original_Seq_Length, Mature_Seq_Length. Mapped UniProtKB_AC. Processed both the main dataset and the initially filtered-out set.

Homology Searching & Structural Context:

Processed previous search results (results_vs_pdb_v2.m8, afdb_found_uniprot_acs_or_upi.csv) to flag proteins with PDB/AFDB hits (Has_Known_Structure).

Performed screen against AFESM "ESM-only" clusters using UniProt IDs.

Running MMseqs2 search against MGnify database on AWS EC2 (Status: Running).

Exploratory Data Analysis: Used Jupyter notebooks (notebooks/) with Python (pandas, matplotlib, seaborn, arcadia-pycolor) to analyze dataset composition, sequence features, functional annotations, and intra-OG sequence divergence.

(Note: Detailed phylogenetic tree inference is being conducted as part of a separate, parallel project.)


> Example:
>
> 1.  Download scripts using `download.ipynb`.
> 2.  Preprocess using `./preprocessing.sh -a data/`
> 3.  Run Snakemake pipeline `snakemake --snakefile Snakefile`
> 4.  Generate figures using `pub/make_figures.ipynb`.

### Compute Specifications

Local: Apple MacBook Pro M3 Max, 36 GB RAM (Used for scripting, USPNet, Metapredict, initial InterProScan runs, data analysis, Git management).

Cloud (AWS EC2):

m6a.24xlarge (96 vCPU, 384 GiB RAM): Used for DeepTMHMM run attempt.

r6i.16xlarge (64 vCPU, 512 GiB RAM): Used for MGnify database indexing and MMseqs2 search.

Operating System (EC2): Amazon Linux 2

Key Software: Python 3.9 (via Conda/Mamba), Pandas, NumPy, Matplotlib, Seaborn, Biopython, MMseqs2 (v17+), InterProScan (v5.73-104.0 via Docker), DeepTMHMM (via Docker), USPNet (local install), Metapredict, OrthoFinder, Docker, AWS CLI, CD-HIT, Git.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please adjust or edit this section as appropriate when you're ready to share your repo.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### VSCode
This template includes recommendations to VSCode users for extensions, particularly the `ruff` linter. These recommendations are stored in `.vscode/extensions.json`. When you open the repository in VSCode, you should see a prompt to install the recommended extensions.

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### `pyproject.toml`
`pyproject.toml` is a configuration file to specify your project's metadata and to set the behavior of other tools such as linters, type checkers etc. You can learn more [here](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)

### Linting
This template automates linting and formatting using GitHub Actions and the `ruff` linter. When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved.
