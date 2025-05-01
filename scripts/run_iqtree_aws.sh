#!/bin/bash

# === AWS EC2 Parallel IQ-TREE Execution Script ===

# --- Configuration ---

# Base working directory on the attached EBS volume
BASE_DIR="/mnt/data" # MODIFY if you used a different mount point

# Directory containing the extracted input alignment files (from trimal_gappyout/)
# Assumes you extracted 'trimal_gappyout_alignments.tar.gz' inside BASE_DIR
INPUT_DIR="${BASE_DIR}/trimal_gappyout/"

# Directory to save the IQ-TREE output files
OUTPUT_DIR="${BASE_DIR}/iqtree_results_gappyout/"

# Directory to save logs (including the parallel job log)
LOG_DIR="${BASE_DIR}/iqtree_logs/"

# Number of parallel IQ-TREE jobs to run simultaneously using GNU Parallel
# *** ADJUST THIS based on your EC2 instance's vCPUs and the -T setting below ***
# Example: Instance with 32 vCPUs, running IQ-TREE with -T 4 => N_JOBS=8 (32/4)
# Example: Instance with 16 vCPUs, running IQ-TREE with -T 4 => N_JOBS=4 (16/4)
N_JOBS=16 # MODIFY AS NEEDED based on instance vCPUs

# Path to IQ-TREE executable (ensure it's installed and in PATH or provide full path)
IQTREE_EXE="iqtree2" # Or "/path/to/iqtree2-X.X.X-Linux/bin/iqtree2"

# IQ-TREE arguments (excluding input -s and output --prefix)
# Updated based on Project Status (Apr 23): MFP, UFBoot, SH-aLRT, 4 threads
# Note: -T 4 allocates 4 threads PER IQ-TREE JOB run by parallel.
IQTREE_ARGS="-m MFP -B 1000 -alrt 1000 -T 4 -redo -quiet"

# Path to GNU Parallel executable (ensure it's installed)
PARALLEL_EXE="parallel"

# Expected input alignment file extension (e.g., .aln, .phy, .fasta)
# *** MODIFY THIS if your files in trimal_gappyout/ have a different extension ***
INPUT_EXT=".aln"

# --- Script Logic ---

echo "Starting IQ-TREE phylogenetic inference on AWS EC2 using GNU Parallel..."
echo "Base directory: ${BASE_DIR}"
echo "Input Alignment directory: ${INPUT_DIR}"
echo "Output IQ-TREE directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo "IQ-TREE Args: ${IQTREE_ARGS}"
echo "Running up to ${N_JOBS} IQ-TREE jobs concurrently (each using ${IQTREE_ARGS} threads)."

# Create output directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Check if IQ-TREE executable is found
if ! command -v "${IQTREE_EXE}" &> /dev/null && [ ! -f "${IQTREE_EXE}" ];
then
    echo "Error: IQ-TREE executable '${IQTREE_EXE}' not found in PATH or as a direct file path." >&2
    exit 1
fi

# Check if GNU Parallel executable is found
if ! command -v "${PARALLEL_EXE}" &> /dev/null
then
    echo "Error: GNU Parallel executable '${PARALLEL_EXE}' not found. Please install it (e.g., 'sudo apt install parallel' or 'sudo yum install parallel')." >&2
    exit 1
fi

# Define the log file for GNU Parallel job tracking
PARALLEL_LOG="${LOG_DIR}/parallel_iqtree_jobs.log"
echo "GNU Parallel job log will be saved to: ${PARALLEL_LOG}"

# Find input files based on INPUT_DIR and INPUT_EXT, count them
input_alignments_list=$(find "${INPUT_DIR}" -maxdepth 1 -name "*${INPUT_EXT}" -print 2>/dev/null)
if [ -z "$input_alignments_list" ]; then
    total_files=0
else
    # Use find again with wc -l for accurate count, handle potential spaces
    total_files=$(find "${INPUT_DIR}" -maxdepth 1 -name "*${INPUT_EXT}" | wc -l)
fi
echo "Found ${total_files} alignment files (*${INPUT_EXT}) to process in ${INPUT_DIR}"

if [[ ${total_files} -eq 0 ]]; then
    echo "Warning: No alignment files (*${INPUT_EXT}) found in ${INPUT_DIR}. Exiting." >&2
    exit 0
fi

# --- Run using GNU Parallel ---
# {} is replaced by the input filename (full path)
# {/.} removes directory path and extension (basename)
# {/} removes the extension only
# {//} removes the directory path only
# --eta provides estimated time remaining
# --joblog logs job status
# Quoting is important for the command string

echo "Launching IQ-TREE jobs via GNU Parallel..."
echo "Check progress with: tail -f ${PARALLEL_LOG}"
echo "Or monitor resource usage with: htop"

find "${INPUT_DIR}" -maxdepth 1 -name "*${INPUT_EXT}" -print0 | \
    "${PARALLEL_EXE}" -0 -j ${N_JOBS} --eta --joblog "${PARALLEL_LOG}" \
    'echo "[$(date)] Starting IQ-TREE for {//}"; '"${IQTREE_EXE}"' -s {} --prefix "'"${OUTPUT_DIR}"'/{/.}" '"${IQTREE_ARGS}"'; echo "[$(date)] Finished IQ-TREE for {//}"'

# Note: Changed prefix slightly to just OG name, IQ-TREE adds suffixes like .iqtree, .log etc.

echo "-------------------------------------"
echo "GNU Parallel IQ-TREE process initiated/completed."
echo "Check the main job log '${PARALLEL_LOG}' for status of each job."
echo "Check individual IQ-TREE log files (e.g., *.log) in ${OUTPUT_DIR}"
echo "IQ-TREE results (.treefile, etc.) will be in: ${OUTPUT_DIR}"
echo "-------------------------------------"

exit 0