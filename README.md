**2025-Asgard-giant-virus-proteomes**

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This repository contains the scripts and analysis notebooks used to construct and analyze a deeply annotated  proteome database for Asgard archaea and Giant Viruses (GV, or NCLDV). The primary motivation is to leverage integrated functional, topological, localization, and structural context predictions, alongside homology searches, to:

Facilitate novel discovery into generalizable "rules" of protein sequence/structure/function relationships.

Identify and prioritize structurally uncharacterized proteins for future study.

Identify orthologous groups relevant to the evolutionary relationships between Asgard archaea and eukaryotes (detailed phylogenetic analysis is part of a parallel project).

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda here. After installing conda and mamba, run the following command to create the primary pipeline run environment used for most analyses. The `envs/dev.yml` file contains the list of dependencies for this environment. Ensure this file is kept up-to-date as new dependencies are added (see instructions below).

# Environment used for data processing, analysis notebooks, MMseqs2 etc.
mamba env create -n asgard_gv_env --file envs/dev.yml
conda activate asgard_gv_env

For specific tools like DeepTMHMM and USPNet, separate environments or setups (e.g., Docker, specific Python versions) may be required. For USPNet, refer to its repository for installation instructions: [ml4bio/USPNet](https://github.com/ml4bio/USPNet).

Install your pre-commit hooks:
```bash
pre-commit install
```
This installs the pre-commit hooks defined in your config (./.pre-commit-config.yaml).

Export your conda environment before sharing:

As your project develops, the number of dependencies in your environment may increase. Whenever you install new dependencies (using either pip install or mamba install), you should update the environment file using the following command.
```bash
conda env export --no-builds > envs/dev.yml
```
`--no-builds` removes build specification from the exported packages to increase portability between different platforms.

## Data

## Input Data:

Source Proteomes: FASTA files for Asgard archaea (311 proteomes), Giant Viruses (446 proteomes), Eukaryotes (63 proteomes), TACK/Euryarchaeota (outgroups). Accessions at data/reference/genome_assembly_list.csv.

Reference Data: interpro_entry.txt (from InterPro), UniProtKB/RefSeq data used implicitly by tools, AFDB_seq_db, PDB_seq_db

Mapping Files: integrated_asgard_gv_ortho_interpro.parquet.

## Intermediate Data:

Filtered FASTA files (e.g., data/input/filtered_fastas.Asgard_all_globular_proteins.fasta).

OrthoFinder output directories (e.g., data/orthofinder_results_Asgard_Orthofinder_Results, data/orthofinder_results/GV_Orthofinder_Results).

InterProScan output directories (TSV, GFF3, e.g., data/interproscan_result_files/Asgard_annotated_nr90.tsv, ).

USPNet intermediate/output directories (eg. USPNet_Processed_Data*/results.csv).

AlphaFold DB search files: afdb_found_uniprot_acs_or_upi.csv

MMseqs2 database indices (e.g., MGnify_DB/).

MMseqs2 search results (results_vs_mgnify.m8).

PDB sequence homology Search Results: results_vs_pdb_v2.m8.

Extracted sequence/ID lists (e.g., unique_virus_names.txt, afesm_esm_only_uniprot_ids.txt).

## Primary Output Data:

proteome_database_v3.5.csv: Latest versions of the main integrated database. 

output_plots/: Directory containing generated figures from analysis notebooks.

## Data Deposition: (TODO: Add Zenodo DOI or other repository links if data is deposited).

Zenodo (tbd)

## Example Workflow / Running the Pipeline

This section outlines a potential workflow for using the scripts in the `scripts/` directory. Users will need to adjust paths, input/output names, and specific parameters according to their data and system setup. It's recommended to consult individual script `--help` messages for detailed options.

**Assumptions:**
*   You have activated the appropriate conda environment (e.g., `asgard_gv_env`).
*   External tools (DIAMOND, InterProScan, MAFFT, FastTree, USPNet, OrthoFinder) are installed and accessible in your PATH, or their paths are provided to the respective wrapper scripts.
*   Input FASTA files are organized appropriately (e.g., one directory per genome for `process_input_faa.py`, or a single directory of FASTA files for others).

---
**Phase 1: Initial Data Preparation & Preprocessing**
---

1.  **Process NCBI FASTA files (if applicable):**
    *   Standardizes headers, extracts bracketed info (e.g., taxonomy).
    *   Script: `scripts/process_input_faa.py`
    *   Example:
        ```bash
        python scripts/process_input_faa.py \
            -i path/to/raw_genome_faa_directories/ \
            -o data/processed_fastas/step1_standardized/ \
            --input_faa_name protein.faa
        ```
    *   Input: Directory containing subdirectories, each with a `protein.faa` (or specified name).
    *   Output: Directory with processed FASTA files, one per input genome, with standardized headers.

2.  **Filter sequences by length:**
    *   Script: `scripts/fasta_length_filter.py`
    *   Example:
        ```bash
        python scripts/fasta_length_filter.py \
            -i data/processed_fastas/step1_standardized/ \
            -o data/processed_fastas/step2_len_filtered/ \
            --min_len 80 \
            --max_len 10000
        ```
    *   Input: Directory of FASTA files (e.g., from step 1).
    *   Output: Directory with FASTA files containing only sequences within the specified length range.

3.  **Filter by predicted disorder (Metapredict):**
    *   Script: `scripts/run_metapredict_filter.py`
    *   Example:
        ```bash
        python scripts/run_metapredict_filter.py \
            -i data/processed_fastas/step2_len_filtered/ \
            -g data/processed_fastas/step3_globular_filtered/ \
            -x data/processed_fastas/step3_skipped_X.fasta \
            -d data/processed_fastas/step3_disordered.fasta \
            -t 0.5
        ```
    *   Input: Directory of FASTA files (e.g., from step 2).
    *   Output:
        *   Directory `step3_globular_filtered/` with FASTA files of predicted globular proteins. This is typically used as input for OrthoFinder.
        *   `step3_skipped_X.fasta` for sequences with 'X' amino acids.
        *   `step3_disordered.fasta` for predicted disordered sequences.

4.  **(Optional) Concatenate FASTA files:**
    *   Concatenate all filtered proteomes (e.g., from `step3_globular_filtered/`) if a single file is needed for tools like InterProScan or USPNet.
    *   Script: `scripts/cat_filter_fastas.py`
    *   Example (concatenating all globular proteins from Step 3):
        ```bash
        python scripts/cat_filter_fastas.py \
            -i data/processed_fastas/step3_globular_filtered/ \
            -c data/processed_fastas/all_globular_proteins_concat.fasta \
            -f data/processed_fastas/all_globular_proteins_dummy_filter.fasta # Filter output can be a dummy if only concat needed
        ```
    *   Output: A single concatenated FASTA file (e.g., `all_globular_proteins_concat.fasta`).

---
**Phase 2: Orthology Analysis (OrthoFinder)**
---
*This is a primary analysis performed after initial proteome filtering.*

5.  **Run OrthoFinder (External Tool):**
    *   Input: Per-genome FASTA files from Step 3 (e.g., `data/processed_fastas/step3_globular_filtered/`).
    *   Example (conceptual):
        ```bash
        orthofinder -f data/processed_fastas/step3_globular_filtered/ -t 16 -o data/orthofinder_results/
        # Or using the Docker image mentioned in the main README:
        # docker run --rm -v /path/to/data:/data davidemms/orthofinder orthofinder -f /data/processed_fastas/step3_globular_filtered/ -t 16 -o /data/orthofinder_results/
        ```
    *   Output: An OrthoFinder results directory (e.g., `data/orthofinder_results/OrthoFinder/Results_*/`).

---
**Phase 3: Functional Annotation (InterProScan)**
---
*This is another primary analysis, run on the comprehensive set of filtered proteins.*

6.  **Run InterProScan:**
    *   Input: A FASTA file containing all proteins to be annotated (e.g., `all_globular_proteins_concat.fasta` from optional Step 4).
    *   Script: `scripts/run_interproscan.sh`
    *   Example:
        ```bash
        bash scripts/run_interproscan.sh \
            -i data/processed_fastas/all_globular_proteins_concat.fasta \
            -d data/interproscan_results/ \
            -f TSV,GFF3 \
            -goterms \
            -pa \
            --cpu 16
        ```
    *   Output: Directory with InterProScan results (TSV, GFF3, XML, etc.).

---
**Phase 4: OrthoFinder Results Analysis & Integration**
---
*Analyzing orthogroups and integrating with metadata. This often uses outputs from Phase 1, 2, and 3.*

7.  **Analyze OrthoFinder results:**
    *   Script: `scripts/orthofinder_analysis.py`
    *   Example:
        ```bash
        python scripts/orthofinder_analysis.py \
            -r data/orthofinder_results/OrthoFinder/Results_*/ \
            -f data/processed_fastas/step3_globular_filtered/ \
            --phylum_map data/reference_tables/genome_to_phylum_map.tsv \
            -p MyAnalysis_OGs \
            --plot_dir output_plots/orthofinder_analysis/ \
            --conservation_threshold 80
        ```
    *   Input: OrthoFinder results directory, directory of FASTA files used for OrthoFinder (for metadata), optional phylum map.
    *   Output: TSV files with OG statistics, lists of conserved OGs, summary text file, and plots.

---
**Phase 5: Other Annotation, ID Mapping & Homology Searches**
---
*These scripts provide additional annotations or mappings and can be run as needed.*

8.  **Run DIAMOND homology search (optional):**
    *   Script: `scripts/run_diamond_loop.sh`
    *   Example (searching individual genome files from Step 3):
        ```bash
        bash scripts/run_diamond_loop.sh \
            -i data/processed_fastas/step3_globular_filtered/ \
            -o data/diamond_results/ \
            -d path/to/your/uniref100.dmnd \
            -t 16 -m 5 -p "*.fasta"
        ```
    *   Output: Directory with DIAMOND search result files.

9.  **Select outgroups from DIAMOND results (for specific phylogenies):**
    *   Script: `scripts/select_outgroups_from_diamond.py`
    *   Example:
        ```bash
        python scripts/select_outgroups_from_diamond.py \
            -i data/diamond_og_vs_outgroups_db_results/ \
            -o data/outgroup_selection/selected_outgroups.csv \
            --pattern "OG*_hits.tsv"
        ```
    *   Output: A CSV file listing selected outgroup sequence IDs.

10. **Predict signal peptides (USPNet):**
    *   Script: `scripts/run_uspnet.sh`
    *   Example (using concatenated FASTA from optional Step 4):
        ```bash
        bash scripts/run_uspnet.sh \
            -i data/processed_fastas/all_globular_proteins_concat.fasta \
            -p USPNet_Intermediate_Data/ \
            -u path/to/cloned/USPNet_scripts/
        ```
    *   Output: CSV file with predictions.

11. **Fetch pLDDT scores from AlphaFold DB:**
    *   Script: `scripts/fetch_plddt.py`
    *   Example:
        ```bash
        python scripts/fetch_plddt.py \
            -i data/reference_tables/my_protein_to_uniprot_map.csv \
            -o data/plddt_scores/protein_plddt_scores.csv \
            --protein_id_col MyProteinID \
            --afdb_id_col CorrespondingUniprotAC_or_UPI
        ```
    *   Output: CSV file with pLDDT scores.

12. **Map UniProtKB ACs to PDB IDs:**
    *   Script: `scripts/uniprot_pdb_search.py`
    *   Example:
        ```bash
        python scripts/uniprot_pdb_search.py \
            -i data/uniprot_ids_for_pdb_search.txt \
            -o data/pdb_mapping/uniprot_to_pdb.tsv
        ```
    *   Output: TSV file mapping UniProtKB ACs to PDB IDs.

13. **Map various protein IDs to UniParc IDs:**
    *   Script: `scripts/uniparc_search.py`
    *   Example:
        ```bash
        python scripts/uniparc_search.py \
            -i data/other_ids_for_uniparc_search.txt \
            -o data/uniparc_mapping/ \
            --output_prefix other_ids
        ```
    *   Output: TSV mapping, not found list, and error log.

---
**Phase 6: Downstream Analyses on Orthogroups**
---
*These steps typically follow orthogroup identification and analysis (Phase 2 & 4).*

14. **Extract sequences for specific OGs:**
    *   Script: `scripts/extract_sequences_by_id.py`
    *   Input: A comprehensive reference FASTA (e.g., `all_globular_proteins_concat.fasta` from optional Step 4) and per-OG ID lists.
    *   Example:
        ```bash
        python scripts/extract_sequences_by_id.py \
            -r data/processed_fastas/all_globular_proteins_concat.fasta \
            -d data/og_id_lists_for_extraction/ \
            -o data/og_sequences_for_alignment/ \
            --hits_suffix .txt \
            --hits_column ProteinID
        ```
    *   Output: Directory with FASTA files, one per OG.

15. **Align sequences within OGs (MAFFT):**
    *   Script: `scripts/run_mafft_parallel.py`
    *   Input: Directory of per-OG FASTA files (e.g., from Step 14).
    *   Example:
        ```bash
        python scripts/run_mafft_parallel.py \
            -i data/og_sequences_for_alignment/ \
            -o data/og_alignments/ \
            -l data/og_alignments_logs/ \
            --input_suffix .fasta \
            --output_suffix .mafft.fa \
            -n 16
        ```
    *   Output: Directory of aligned FASTA files.

16. **Build phylogenetic trees from alignments (FastTree):**
    *   Script: `scripts/run_fasttree_parallel.py`
    *   Input: Directory of aligned FASTA files (e.g., from Step 15).
    *   Example:
        ```bash
        python scripts/run_fasttree_parallel.py \
            -i data/og_alignments/ \
            -o data/og_trees/ \
            -l data/og_trees_logs/ \
            --input_suffix .mafft.fa \
            --output_suffix .tree \
            -n 16
        ```
    *   Output: Directory of tree files (Newick format).

17. **Perform Hill diversity analysis (example):**
    *   Script: `scripts/hill_diversity_analysis.py`
    *   Example:
        ```bash
        python scripts/hill_diversity_analysis.py \
            -i data/analysis_inputs/og_domain_counts.csv \
            -o output_plots/hill_diversity_domains.png
        ```
    *   Output: PNG plot of Hill diversity.

---

This workflow provides a comprehensive guide. Remember to adapt file paths, names, and specific tool parameters based on your actual data and analytical goals. For detailed options of each script, run `python scripts/script_name.py --help`.

## Overview

### Description of the folder structure


### Methods

Data Acquisition & Preparation: Downloaded source proteomes (Asgard: 311, GV: 451) from NCBI. Filtered proteomes for length (80-1000 aa) and predicted globularity (Metapredict). Generated custom FASTA headers.

Orthology Inference: Ran OrthoFinder separately on the filtered Asgard and GV proteomes to assign proteins to Orthogroups (OGs). A Docker image for OrthoFinder can be found here: [davidemms/OrthoFinder-Dockerfile](https://github.com/davidemms/OrthoFinder-Dockerfile).

Functional Annotation (InterProScan): Ran InterProScan (v5.73-104.0 via Docker) iteratively. Used a curated application list (Pfam, CDD, Gene3D, SUPERFAMILY, SMART, ProSitePatterns, ProSiteProfiles) and requested GO terms/pathways. Managed memory via modified interproscan.sh.

Custom Functional Categorization: Developed and applied a rule-based Python script (add_specific_category_IPR_v10.py) using specific IPR IDs, keywords (in IPR names/source annotations), and IPR type fallback to assign Specific_Functional_Category and Category_Trigger.

Taxonomy Refinement: Fetched NCBI TaxIDs via Entrez; parsed Asgard Phylum, Virus Name; assigned Virus Family via custom keyword matching.

Sequence Feature Prediction:

Ran USPNet locally to predict signal peptides and infer likely localization (Signal_Peptide_USPNet, SP_Cleavage_Site_USPNet).

Ran Metapredict locally via add_disorder.py script to calculate Percent_Disorder.

Database Integration: Iteratively merged all annotations and metadata into the proteome_database CSV using the database_assembly.ipynb notebook. 

Homology Searching & Structural Context:

Processed previous search results (results_vs_pdb_v2.m8, afdb_found_uniprot_acs_or_upi.csv) to flag proteins with PDB/AFDB hits (Has_Known_Structure).

Performed screen against AFESM "ESM-only" clusters using UniProt IDs. The AFESM "ESM-only" clusters were derived from database associated with "Metagenomic-scale analysis of the predicted protein structure universe; https://doi.org/10.1101/2025.04.23.650224". First, we ran all Structurally Dark proteins on an MMseqs2 search against the MGnify database, to identify MGnify clusters these proteins fall into. Then, we searched these against the AFESM "ESM-only" clusters, to determine if any of the Structurall Dark proteins were in the ESMAtlas.

Exploratory Data Analysis: Used Jupyter notebooks (notebooks/) with Python (pandas, matplotlib, seaborn, arcadia-pycolor) to analyze dataset composition, sequence features, functional annotations, and intra-OG sequence divergence.

### Compute Specifications

Local: Apple MacBook Pro M3 Max, 36 GB RAM (Used for scripting, USPNet, Metapredict, InterProScan runs, data analysis, Git management).

Cloud (AWS EC2):

r6i.16xlarge (64 vCPU, 512 GiB RAM): Used for MGnify database indexing and MMseqs2 search.

Operating System (EC2): Amazon Linux 2

Key Software: Python 3.9 (via Conda/Mamba), Pandas, NumPy, Matplotlib, Seaborn, Biopython, MMseqs2 (v17+), InterProScan (v5.73-104.0 via Docker), USPNet (local install), Metapredict, OrthoFinder, Docker, AWS CLI.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please adjust or edit this section as appropriate when you're ready to share your repo.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### `pyproject.toml`
`pyproject.toml` is a configuration file to specify your project's metadata and to set the behavior of other tools such as linters, type checkers etc. You can learn more [here](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)

### Linting
This template automates linting and formatting using GitHub Actions and the `ruff` linter. When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved.
