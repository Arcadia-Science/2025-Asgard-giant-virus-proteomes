**2025-Asgard-giant-virus-proteomes**

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This repository contains the scripts and analysis notebooks used to construct and analyze a deeply annotated  proteome database for Asgard archaea and Giant Viruses (GV, or NCLDV). The primary motivation is to leverage integrated functional, topological, localization, and structural context predictions, alongside homology searches, to:

Facilitate novel discovery into generalizable "rules" of protein sequence/structure/function relationships.

Identify and prioritize structurally uncharacterized proteins for future study.

Identify orthologous groups relevant to the evolutionary relationships between Asgard archaea and eukaryotes (detailed phylogenetic analysis is part of a parallel project).

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda here. After installing conda and mamba, run the following command to create the primary pipeline run environment used for most analyses. Note: Specific tools like DeepTMHMM and USPNet required separate environments/setups (Docker, specific Python versions).

# Environment used for data processing, analysis notebooks, MMseqs2 etc.

mamba env create -n asgard_gv_env --file envs/dev.yml
conda activate asgard_gv_env

The `envs/dev.yml` file contains the list of dependencies for the primary Conda environment. Ensure this file is kept up-to-date as new dependencies are added (see below).

For specific tools like DeepTMHMM and USPNet, separate environments or setups (e.g., Docker, specific Python versions) may be required. For USPNet, refer to its repository for installation instructions: [ml4bio/USPNet](https://github.com/ml4bio/USPNet).

Install your pre-commit hooks:

pre-commit install

This installs the pre-commit hooks defined in your config (./.pre-commit-config.yaml).

Export your conda environment before sharing:

As your project develops, the number of dependencies in your environment may increase. Whenever you install new dependencies (using either pip install or mamba install), you should update the environment file using the following command.

conda env export --no-builds > envs/dev.yml
```--no-builds` removes build specification from the exported packages to increase portability between different platforms.

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
*   You have activated an appropriate conda environment.
*   External tools (DIAMOND, InterProScan, MAFFT, FastTree (or VeryFastTree), USPNet, OrthoFinder) are installed and accessible in your PATH, or their paths are provided to the respective wrapper scripts.
*   Input FASTA files are organized appropriately (e.g., one directory per genome for `process_input_faa.py`, or a single directory of FASTA files for others).

---
**Phase 1: Initial Data Preparation & Preprocessing**
---

1.  **Process NCBI FASTA files (if applicable):**
    *   Standardizes headers, extracts bracketed info (e.g., taxonomy).
    *   Script: `scripts/process_input_faa.py`
    *   Example:
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
        python scripts/run_metapredict_filter.py \
            -i data/processed_fastas/step2_len_filtered/ \
            -g data/processed_fastas/step3_globular_filtered/ \
            -x data/processed_fastas/step3_skipped_X.fasta \
            -d data/processed_fastas/step3_disordered.fasta \
            -t 0.5 
        ```
    *   Input: Directory of FASTA files (e.g., from step 2).
    *   Output: 
        *   Directory `step3_globular_filtered/` with FASTA files of predicted globular proteins.
        *   `step3_skipped_X.fasta` for sequences with 'X' amino acids.
        *   `step3_disordered.fasta` for predicted disordered sequences.
    *   *Use `data/processed_fastas/step3_globular_filtered/` for subsequent steps if selecting for globular proteins.*

4.  **Concatenate and/or filter by keywords (optional, example for hypotheticals):**
    *   This script can be used to create specific subsets, e.g., all hypothetical proteins for InterProScan.
    *   Script: `scripts/cat_filter_fastas.py`
    *   Example (extracting hypotheticals from globular proteins):
        python scripts/cat_filter_fastas.py \
            -i data/processed_fastas/step3_globular_filtered/ \
            -c data/processed_fastas/all_globular_concat.fasta \
            -f data/processed_fastas/hypothetical_globular_subset.fasta \
            --keywords hypothetical uncharacterized "unknown function" \
            --name_index 4 # Adjust if your standardized headers have annotation at a different field
        ```
    *   Input: Directory of FASTA files.
    *   Output: A concatenated FASTA and a filtered FASTA subset.

---
**Phase 2: Orthology Analysis (using external OrthoFinder)**
---

5. **Run OrthoFinder:**
    *   This step uses the OrthoFinder tool, which is external to this script suite.
    *   Input to OrthoFinder would typically be the processed, filtered FASTA files (one per species/genome), e.g., from `data/processed_fastas/step3_globular_filtered/` or a non-redundant set like `all_globular_nr90.fasta` if analyzing a single combined proteome.
    *   Example:
        orthofinder -f data/processed_fastas/step3_globular_filtered/ -t 16 -o data/orthofinder_results/
        # Or using the Docker image mentioned in the main README
        # docker run --rm -v /path/to/data:/data davidemms/orthofinder orthofinder -f /data/processed_fastas/step3_globular_filtered/ -t 16 -o /data/orthofinder_results/
        ```

6. **Analyze OrthoFinder results:**
    *   Script: `scripts/orthofinder_analysis.py`
    *   Example:
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
**Phase 3: Homology Search & Outgroup Selection**
---

7.  **Run DIAMOND homology search:**
    *   Search your processed proteomes against a reference database (e.g., UniRef100).
    *   Script: `scripts/run_diamond_loop.sh`
    *   Example:
        bash scripts/run_diamond_loop.sh \
            -i data/processed_fastas/step3_globular_filtered/ \
            -o data/diamond_results/ \
            -d path/to/your/uniref100.dmnd \
            -t 16 \
            -m 5 \
            -p "*.fasta" # Adjust pattern if your files from step3 have a different extension
        ```
    *   Input: Directory of FASTA files (one per genome/proteome).
    *   Output: Directory with DIAMOND search result files (TSV format), one per input FASTA.

8.  **Select outgroups from DIAMOND results (if needed for phylogenies):**
    *   This script processes DIAMOND results (e.g., against TACK/Euryarchaeota) to pick outgroups for OGs.
    *   *This step often requires DIAMOND searches of OG sequences against a database containing potential outgroups, which is not explicitly generated by `run_diamond_loop.sh` in the example above. Assume you have such DIAMOND results per OG.*
    *   Script: `scripts/select_outgroups_from_diamond.py`
    *   Example:
        python scripts/select_outgroups_from_diamond.py \
            -i data/diamond_og_vs_outgroups_db_results/ \
            -o data/outgroup_selection/selected_outgroups.csv \
            --pattern "OG*_hits.tsv" \
            --evalue 1e-10 \
            --min_cov 0.5 \
            --max_outgroups 3
        ```
    *   Input: Directory of DIAMOND TSV files, where each file contains hits for one OG.
    *   Output: A CSV file listing selected outgroup sequence IDs for each OG.

---
**Phase 4: Functional Annotation & ID Mapping**
---

9.  **Run InterProScan for functional annotation:**
    *   Typically run on a comprehensive, non-redundant set of proteins.
    *   Script: `scripts/run_interproscan.sh`
    *   Example (using the nr90 set from CD-HIT):
        # Ensure interproscan.properties is configured correctly, especially bin.directory
        # and java.command if not in PATH.
        # The script passes arguments directly to InterProScan's interproscan.sh
        bash scripts/run_interproscan.sh \
            -i data/processed_fastas/all_globular_nr90.fasta \
            -d data/interproscan_results/ \
            -f TSV,GFF3 \
            -goterms \
            -pa \
            --cpu 16 
        # Add other InterProScan options as needed (e.g., -appl Pfam,CDD)
        ```
    *   Input: A FASTA file of proteins to annotate.
    *   Output: Directory with InterProScan results (TSV, GFF3, XML, etc.).

10.  **Predict signal peptides (USPNet):**
    *   Script: `scripts/run_uspnet.sh`
    *   Example (using a concatenated FASTA of all proteins of interest):
        # Ensure all_proteins_for_uspnet.fasta is prepared
        bash scripts/run_uspnet.sh \
            -i data/processed_fastas/all_globular_concat.fasta \
            -p USPNet_Intermediate_Data/ \
            -u path/to/cloned/USPNet_scripts/ # If USPNet not installed in PATH
        ```
    *   Input: A single FASTA file.
    *   Output: CSV file with predictions inside the processed data directory (e.g., `USPNet_Intermediate_Data/results.csv`).

11. **Fetch pLDDT scores from AlphaFold DB:**
    *   Requires a CSV input mapping your protein IDs to UniProt ACs or UPIs.
    *   Script: `scripts/fetch_plddt.py`
    *   Example:
        python scripts/fetch_plddt.py \
            -i data/reference_tables/my_protein_to_uniprot_map.csv \
            -o data/plddt_scores/protein_plddt_scores.csv \
            --protein_id_col MyProteinID \
            --afdb_id_col CorrespondingUniprotAC_or_UPI \
            --delay 0.2
        ```
    *   Input: CSV file with protein IDs and AFDB queryable IDs.
    *   Output: CSV file with protein IDs and their average pLDDT scores.

12. **Map UniProtKB ACs to PDB IDs:**
    *   Script: `scripts/uniprot_pdb_search.py`
    *   Example:
        # Create a file with one UniProtKB AC per line
        # e.g., data/uniprot_ids_for_pdb_search.txt
        python scripts/uniprot_pdb_search.py \
            -i data/uniprot_ids_for_pdb_search.txt \
            -o data/pdb_mapping/uniprot_to_pdb.tsv \
            --delay 0.2
        ```
    *   Input: Text file with UniProtKB Accessions.
    *   Output: TSV file mapping UniProtKB ACs to PDB IDs.

13. **Map various protein IDs to UniParc IDs:**
    *   Useful for IDs not found in UniProtKB.
    *   Script: `scripts/uniparc_search.py`
    *   Example:
        # Create a file with one protein ID (e.g., RefSeq) per line
        # e.g., data/other_ids_for_uniparc_search.txt
        python scripts/uniparc_search.py \
            -i data/other_ids_for_uniparc_search.txt \
            -o data/uniparc_mapping/ \
            --output_prefix other_ids 
        ```
    *   Input: Text file with protein IDs.
    *   Output: Directory with a TSV mapping, a list of not found IDs, and an error log.

---
**Phase 5: Downstream Analyses on Orthogroups**
---

14. **Extract sequences for specific OGs:**
    *   Script: `scripts/extract_sequences_by_id.py`
    *   This script needs a "reference FASTA" containing all sequences and "hit files" listing IDs per OG.
    *   The reference FASTA could be `data/processed_fastas/all_globular_concat.fasta`.
    *   Hit files (e.g., `OG00001.txt` containing protein IDs for that OG) would need to be generated, perhaps from `orthofinder_analysis.py` outputs or OrthoFinder's own files.
    *   Example (conceptual, assuming hit files are prepared):
        python scripts/extract_sequences_by_id.py \
            -r data/processed_fastas/all_globular_concat.fasta \
            -d data/og_id_lists_for_extraction/ \
            -o data/og_sequences_for_alignment/ \
            --hits_suffix .txt \
            --hits_column ProteinID 
        ```
    *   Input: Reference FASTA, directory of ID list files.
    *   Output: Directory with FASTA files, one per input ID list (per OG).

15. **Align sequences within OGs (MAFFT):**
    *   Script: `scripts/run_mafft_parallel.py`
    *   Example:
        python scripts/run_mafft_parallel.py \
            -i data/og_sequences_for_alignment/ \
            -o data/og_alignments/ \
            -l data/og_alignments_logs/ \
            --input_suffix .fasta \
            --output_suffix .mafft.fa \
            -n 16 
        ```
    *   Input: Directory of per-OG FASTA files (e.g., from step 15).
    *   Output: Directory of aligned FASTA files (e.g., `OG123.mafft.fa`). Per-job logs.

16. **Build phylogenetic trees from alignments (FastTree):**
    *   Script: `scripts/run_fasttree_parallel.py`
    *   Example:
        python scripts/run_fasttree_parallel.py \
            -i data/og_alignments/ \
            -o data/og_trees/ \
            -l data/og_trees_logs/ \
            --input_suffix .mafft.fa \
            --output_suffix .tree \
            -n 16
        ```
    *   Input: Directory of aligned FASTA files (e.g., from step 16).
    *   Output: Directory of tree files (Newick format). Per-job logs.

17. **Perform Hill diversity analysis (example):**
    *   Requires an input CSV where rows are OGs and columns are features (e.g., domain counts per OG). This CSV would need to be generated from other data (e.g., InterProScan results merged with OG assignments).
    *   Script: `scripts/hill_diversity_analysis.py`
    *   Example (assuming `og_domain_counts.csv` is prepared):
        python scripts/hill_diversity_analysis.py \
            -i data/analysis_inputs/og_domain_counts.csv \
            -o output_plots/hill_diversity_domains.png
        ```
    *   Input: CSV file with orthogroups and feature counts/proportions.
    *   Output: PNG plot of Hill diversity.

---

This workflow provides a comprehensive guide. Remember to adapt file paths, names, and specific tool parameters based on your actual data and analytical goals. For detailed options of each script, run `python scripts/script_name.py --help`.

## Overview

### Description of the folder structure


### Methods

Data Acquisition & Preparation: Downloaded source proteomes (Asgard: 311, GV: 451) from NCBI. Filtered proteomes for length (80-1000 aa) and predicted globularity (Metapredict). Generated custom FASTA headers.

Orthology Inference: Ran OrthoFinder separately on the filtered Asgard and GV proteomes to assign proteins to Orthogroups (OGs). A Docker image for OrthoFinder can be found here: [davidemms/OrthoFinder-Dockerfile](https://github.com/davidemms/OrthoFinder-Dockerfile).

Functional Annotation (InterProScan):

Ran InterProScan (v5.73-104.0 via Docker) iteratively. Initially on nr90 subsets (created via cd-hit -c 0.90), then on proteins initially lacking hits to maximize coverage. Used a curated application list (Pfam, CDD, Gene3D, SUPERFAMILY, SMART, ProSitePatterns, ProSiteProfiles) and requested GO terms/pathways. Managed memory via modified interproscan.sh.

Custom Functional Categorization: Developed and applied a rule-based Python script (add_specific_category_IPR_v10.py) using specific IPR IDs, keywords (in IPR names/source annotations), and IPR type fallback to assign Specific_Functional_Category and Category_Trigger.

Taxonomy Refinement: Fetched NCBI TaxIDs via Entrez; parsed Asgard Phylum, Virus Name; assigned Virus Family via custom keyword matching.

Sequence Feature Prediction:

Ran USPNet locally to predict signal peptides and infer likely localization (Signal_Peptide_USPNet, SP_Cleavage_Site_USPNet).

Ran Metapredict locally via add_disorder.py script to calculate Percent_Disorder to the main database.

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

Key Software: Python 3.9 (via Conda/Mamba), Pandas, NumPy, Matplotlib, Seaborn, Biopython, MMseqs2 (v17+), InterProScan (v5.73-104.0 via Docker), USPNet (local install), Metapredict, OrthoFinder, Docker, AWS CLI, CD-HIT, Git.

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
