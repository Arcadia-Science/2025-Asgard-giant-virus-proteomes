**2025-Asgard-giant-virus-proteomes**

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This repository contains the scripts and analysis notebooks used to construct and analyze a deeply annotated  proteome database for Asgard archaea and Giant Viruses (GV, or NCLDV). The primary motivation is to leverage integrated functional, topological, localization, and structural context predictions, alongside homology searches, to:

- Facilitate novel discovery into generalizable "rules" of protein sequence/structure/function relationships.
- Identify and prioritize structurally uncharacterized proteins for future study.
- Identify orthologous groups relevant to the evolutionary relationships between Asgard archaea and eukaryotes (detailed phylogenetic analysis is part of a parallel project).

### Data Availability

The large data files required to run this workflow (e.g., raw proteomes, DIAMOND databases, and large intermediate files) are not stored in this Git repository. They are permanently archived and publicly available on Zenodo.
**To run the analysis, please first download the data from Zenodo [https://doi.org/10.5281/zenodo.15933361] and place the contents into the corresponding subdirectories within the `data/` folder in this repository.**

After extraction, the `data/` directory should contain the following subdirectories (the Git repository already includes these as empty/partial placeholders):

```
data/
├── input/                 # Raw and primary input data files (e.g., raw proteomes)
├── intermediate/          # Files generated during intermediate processing steps
├── output/                # Final data outputs (analysis results, summary tables)
├── reference/             # Reference files (taxonomies, domain information)
└── orthofinder_results/   # OrthoFinder outputs
    ├── Asgard_OrthoFinder_Results/
    └── GV_Orthofinder_Results/
```

After downloading, verify file integrity against the checksums published on the
Zenodo record page (each file lists an MD5 hash). For example:

```bash
md5sum -c <(echo "<md5_from_zenodo>  <downloaded_file>")
```


### Project Structure

This repository is organized into the following directories:

-   **`/data`**: Contains all data files, organized by status.
    -   `input/`: Raw and primary input data files.
    -   `intermediate/`: Files generated during intermediate processing steps.
    -   `output/`: Final data outputs, such as analysis results and summary tables.
    -   `reference/`: Reference files like taxonomies or domain information.
-   **`/envs`**: Contains the Conda environment file (`dev.yml`) for reproducing the computational environment.
-   **`/notebooks`**: Jupyter notebooks for exploratory data analysis, figure generation, and summarizing results. The `_refactored` notebooks are the canonical versions used for the publication; the originals are retained for provenance.
    - `database_assembly_refactored.ipynb` (canonical; original: `database_assembly.ipynb`): Notebook for assembling and analyzing the database.
    - `Figures_DB_Pub_v2_refactored.ipynb` (canonical; original: `Figures_DB_Pub_v2.ipynb`): Notebook for generating the figures for the publication.
-   **`/scripts`**: Contains the core analysis scripts and workflow wrappers essential for the scientific findings of this paper.
    -   `utils/`: Contains helper scripts for generic data pre-processing tasks, such as filtering FASTA files and formatting headers.
      
## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda here. After installing conda and mamba, run the following command to create the primary pipeline run environment used for most analyses. Note: Specific tools like DeepTMHMM and USPNet required separate environments/setups (Docker, specific Python versions).

```bash
# Environment used for data processing, analysis notebooks, MMseqs2 etc.
mamba env create -n asgard_gv_env --file envs/dev.yml
conda activate asgard_gv_env
```

The `envs/dev.yml` file contains the list of dependencies for the primary Conda environment. Ensure this file is kept up-to-date as new dependencies are added (see below).

For USPNet, refer to its repository for installation instructions: [ml4bio/USPNet](https://github.com/ml4bio/USPNet).

## Data

### Reference Data 

- `interpro_entry.txt` (https://www.ebi.ac.uk/interpro/download/)
- `genome_assembly_taxonomy_list.csv` (generated from the downloaded assemblies)
- `eukaryotic_species_taxonomy.csv` (generated from the downloaded eukaryotic proteomes)
- `afdb_protein_confidence_scores.csv` (generated using scripts/fetch_plddt.py)

### Input Data

Source Proteomes: 
- FASTA files for Asgard archaea (311 proteomes)
- Giant Viruses (446 proteomes)
- Eukaryotes (63 proteomes)

Accessions at `data/reference/genome_assembly_list.csv` and `data/reference/eukaryotic_species_taxonomy.csv`

### Intermediate Data

Many of these files are large, and are stored in a Zenodo repository (DOI: 10.5281/zenodo.15933361). They include:

- Filtered FASTA files (e.g., `data/input/filtered_fastas.Asgard_all_globular_proteins.fasta`).
- OrthoFinder output files.
- InterProScan output mapping (`integrated_asgard_gv_ortho_interpro.parquet`).
- USPNet intermediate/output directories (eg. `USPNet_Processed_Data*/results.csv`).
- MMseqs2 search results: `results_vs_mgnify_ESMA.m8`.
- PDB sequence homology Search Results: `results_vs_pdb_v2.m8`.
- AFDB mapping results: `afdb_found_uniprot_acs_or_upi.csv`
- AFESM ESM-only clusters: `afesm_esm_only_uniprot_ids.txt`
- AFESM MGnify clusters: `afesm_mgnify_uniprot_ids.txt`
- Extracted sequence/ID lists (e.g., `unique_virus_names.txt`, `afesm_esm_only_uniprot_ids.txt`).

### Primary Output Data

- `proteome_database_v3.5.csv`: Latest versions of the main integrated database. 
- `output_plots/`: Directory containing generated figures from analysis notebooks.

## Example Workflow / Running the Pipeline

This section outlines a potential workflow for using the scripts in the `scripts/` directory. Users will need to adjust paths, input/output names, and specific parameters according to their data and system setup. It's recommended to consult individual script `--help` messages for detailed options.

**Assumptions:**
*   You have created and activated the conda environment from `envs/dev.yml`, which provides DIAMOND, CD-HIT, MAFFT, and VeryFastTree (the `run_fasttree_parallel.py` wrapper defaults to the `VeryFastTree` executable). The `veryfasttree` package installs its executable as `VeryFastTree`; pass `--fasttree_exe fasttree` if you use the original FastTree instead.
*   Tools not bundled in the conda environment (InterProScan via Docker, USPNet via a local install, and OrthoFinder via Docker) are installed separately and accessible in your PATH, or their paths are provided to the respective wrapper scripts.
*   Input FASTA files are organized appropriately (e.g., one directory per genome for `process_input_faa.py`, or a single directory of FASTA files for others).

### Phase 1: Initial Data Preparation & Preprocessing

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
    *   Script: `scripts/utils/fasta_length_filter.py`
    *   Example:
        ```bash
        python scripts/utils/fasta_length_filter.py \
            -i data/processed_fastas/step1_standardized/ \
            -o data/processed_fastas/step2_len_filtered/ \
            --min_len 80 \
            --max_len 10000
        ```
    *   Input: Directory of FASTA files (e.g., from step 1).
    *   Output: Directory with FASTA files containing only sequences within the specified length range.

3.  **Filter by predicted disorder (Metapredict v3):**
    *   Script: `scripts/run_metapredict_filter.py`
    *   Example:
        ```bash
        python scripts/run_metapredict_filter.py \
            -i data/processed_fastas/step2_len_filtered/ \
            -o data/processed_fastas/step3_globular_filtered/
        ```
    *   Input: Directory of FASTA files (e.g., from step 2).
    *   Output: 
        *   Directory `step3_globular_filtered/` with FASTA files of predicted globular proteins.
        *   `skipped_X.fasta` for sequences with 'X' amino acids.
        *   `sdisordered.fasta` for predicted disordered sequences.
    *   Use `data/processed_fastas/step3_globular_filtered/` for subsequent steps if selecting for globular proteins.

4.  **Concatenate**
    *   This script can be used to create specific subsets, e.g., all hypothetical proteins for InterProScan.
    *   Script: `scripts/utils/cat_filter_fastas.py`
    *   Example (extracting hypotheticals from globular proteins):
        ```bash
        python scripts/utils/cat_filter_fastas.py \
            -i data/processed_fastas/step3_globular_filtered/ \
            -c data/processed_fastas/all_globular_concat.fasta \
            -f data/processed_fastas/all_globular_hypotheticals.fasta
        ```
    *   Input: Directory of FASTA files.
    *   Output: A concatenated FASTA .

### Phase 2: Orthology Analysis (using OrthoFinder docker image)

5. **Run OrthoFinder:**
    *   This step uses the OrthoFinder tool, which is external to this script suite.
    *   Input to OrthoFinder would typically be the processed, filtered FASTA files (one per species/genome), e.g., from `data/processed_fastas/step3_globular_filtered/` or a non-redundant set like `all_globular_nr90.fasta` if analyzing a single combined proteome.
    *   Example:
        ```bash
        docker run --rm -v /path/to/data:/data davidemms/orthofinder \
            orthofinder \
            -f /data/processed_fastas/step3_globular_filtered/ \
            -t 16 \
            -o /data/orthofinder_results/
        ```

6. **Analyze OrthoFinder results:**
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

### Phase 3: Eukaryotic Homology Search

7.  **Run DIAMOND homology search:**
    *   Search your processed proteomes against a reference database (e.g., UniRef100).
    *   Script: `scripts/run_diamond_loop.sh`
    *   Example:
        ```bash
        bash scripts/run_diamond_loop.sh \
            -i data/processed_fastas/step3_globular_filtered/ \
            -o data/diamond_results/ \
            -d path/to/your/db.dmnd \
            -t 16 \
            -m 5 \
            -p "*.fasta"
        ```
    *   Input: Directory of FASTA files (one per genome/proteome).
    *   Output: Directory with DIAMOND search result files (TSV format), one per input FASTA.

### Phase 4: Functional Annotation & ID Mapping

8.  **Run InterProScan for functional annotation:**
    *   Ensure Interproscan and all necessary files are properly installed locally
    *   Script: `scripts/run_interproscan.sh`
    *   Example :
        ```bash
        # The script passes arguments directly to InterProScan's interproscan.sh
        bash scripts/run_interproscan.sh \
            -i data/processed_fastas/Asgard_all_globular_proteins.fasta \
            -d data/interproscan_results/ \
            -f TSV,GFF3 \
            -goterms \
            -pa \
            --cpu 16
        # Add other InterProScan options as needed (e.g., -appl Pfam,CDD)
        ```
    *   Input: A FASTA file of proteins to annotate.
    *   Output: Directory with InterProScan results (TSV, GFF3, XML, etc.).

9.  **Predict signal peptides (USPNet):**
    *   Script: `scripts/run_uspnet.sh`
    *   Example (using a concatenated FASTA of all proteins of interest):
        ```bash
        # Ensure all_proteins_for_uspnet.fasta is prepared
        bash scripts/run_uspnet.sh \
            -i data/processed_fastas/all_globular_concat.fasta \
            -p USPNet_Intermediate_Data/ \
            -u path/to/cloned/USPNet_scripts/ # If USPNet not installed in PATH
        ```
    *   Input: A single FASTA file.
    *   Output: CSV file with predictions inside the processed data directory (e.g., `USPNet_Intermediate_Data/results.csv`).

10. **Fetch pLDDT scores from AlphaFold DB:**
    *   Requires a CSV input mapping your protein IDs to UniProt ACs or UPIs.
    *   Script: `scripts/fetch_plddt.py`
    *   Example:
        ```bash
        python scripts/fetch_plddt.py \
            -i data/reference_tables/my_protein_to_uniprot_map.csv \
            -o data/plddt_scores/protein_plddt_scores.csv \
            --protein_id_col MyProteinID \
            --afdb_id_col CorrespondingUniprotAC_or_UPI \
            --delay 0.2
        ```
    *   Input: CSV file with protein IDs and AFDB queryable IDs.
    *   Output: CSV file with protein IDs and their average pLDDT scores.

11. **Map UniProtKB ACs to PDB IDs:**
    *   Script: `scripts/uniprot_pdb_search.py`
    *   Example:
        ```bash
        # Create a file with one UniProtKB AC per line
        # e.g., data/uniprot_ids_for_pdb_search.txt
        python scripts/uniprot_pdb_search.py \
            -i data/uniprot_ids_for_pdb_search.txt \
            -o data/pdb_mapping/uniprot_to_pdb.tsv \
            --delay 0.2
        ```
    *   Input: Text file with UniProtKB Accessions.
    *   Output: TSV file mapping UniProtKB ACs to PDB IDs.

### Phase 5: Downstream Analyses on Orthogroups

13. **Extract sequences for specific OGs:**
    *   Script: `scripts/extract_sequences_by_id.py`
    *   This script needs a "reference FASTA" containing all sequences and "hit files" listing IDs per OG.
    *   The reference FASTA could be `data/processed_fastas/all_globular_concat.fasta`.
    *   Hit files (e.g., `OG00001.txt` containing protein IDs for that OG) would need to be generated, perhaps from `orthofinder_analysis.py` outputs or OrthoFinder's own files.
    *   Example (conceptual, assuming hit files are prepared):
        ```bash
        python scripts/extract_sequences_by_id.py \
            -r data/processed_fastas/all_globular_concat.fasta \
            -d data/og_id_lists_for_extraction/ \
            -o data/og_sequences_for_alignment/ \
            --hits_suffix .txt \
            --hits_column ProteinID
        ```
    *   Input: Reference FASTA, directory of ID list files.
    *   Output: Directory with FASTA files, one per input ID list (per OG).

14. **Align sequences within OGs (MAFFT):**
    *   Script: `scripts/run_mafft_parallel.py`
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
    *   Input: Directory of per-OG FASTA files (e.g., from step 15).
    *   Output: Directory of aligned FASTA files (e.g., `OG123.mafft.fa`). Per-job logs.

15. **Build phylogenetic trees from alignments (FastTree):**
    *   Script: `scripts/run_fasttree_parallel.py`
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
    *   Input: Directory of aligned FASTA files (e.g., from step 16).
    *   Output: Directory of tree files (Newick format). Per-job logs.

16. **Perform Hill diversity analysis (example):**
    *   Requires an input CSV where rows are OGs and columns are features (e.g., domain counts per OG). This CSV would need to be generated from other data (e.g., InterProScan results merged with OG assignments).
    *   Script: `scripts/hill_diversity_analysis.py`
    *   Example (assuming `og_domain_counts.csv` is prepared):
        ```bash
        python scripts/hill_diversity_analysis.py \
            -i data/analysis_inputs/og_domain_counts.csv \
            -o output_plots/hill_diversity_domains.png
        ```
    *   Input: CSV file with orthogroups and feature counts/proportions.
    *   Output: PNG plot of Hill diversity.

This workflow provides a comprehensive guide. Remember to adapt file paths, names, and specific tool parameters based on your actual data and analytical goals. For detailed options of each script, run `python scripts/script_name.py --help`.

    # Detailed methods of how we used these scripts and analyzed the outputs are available in the Methods section of  the pub (https://doi.org/10.57844/arcadia-prc5-56p7) 

#### Database Integration

All annotations and metadata were iteratively merged into a single `proteome_database` CSV file using the `database_assembly_refactored.ipynb` Jupyter notebook.


### Compute Specifications

Apple MacBook Pro M3 Max, 36 GB RAM (Used for USPNet, Metapredict, InterProScan runs, data analysis).

AWS EC2 (used for MGnify database indexing and MMseqs2 search): r6i.16xlarge (64 vCPU, 512 GiB RAM) with Amazon Linux 2.

Required Software: Python 3.9 (via Conda/Mamba), MMseqs2 (v17+), InterProScan (v5.73-104.0 via Docker), USPNet (local install), Metapredict, OrthoFinder, CD-HIT.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
