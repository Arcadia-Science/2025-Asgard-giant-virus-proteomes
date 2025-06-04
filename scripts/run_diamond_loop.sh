#!/bin/bash

# Script to run DIAMOND blastp for each target OG listed in a file
# against a specified DIAMOND database (e.g., eukaryotic proteomes).

# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# Prevent errors in a pipeline from being masked.
set -euo pipefail

# --- Configuration ---
# REVIEW: Consider using command-line arguments (e.g., getopts) for more flexibility
#         instead of hardcoding paths here, especially if used in different contexts.
TARGET_OG_LIST="analysis_step_3_5_results/target_ogs_for_phylogeny.txt" # File containing one OG ID per line
INPUT_FASTA_DIR="analysis_step_3_5_results/per_og_fastas"              # Directory containing input FASTA files named like ${OG_ID}_phylo_input.fasta
EUK_DB_NAME="euk63_proteomes"                                          # Base name of the DIAMOND DB (e.g., 'euk63_proteomes' for 'euk63_proteomes.dmnd')
OUTPUT_DIR="analysis_step_3_5_results/diamond_search_results"          # Directory to store output TSV files
THREADS=10                                                             # Number of threads for each DIAMOND job
MAX_TARGETS=10                                                         # Max target sequences per query sequence
EVALUE=1e-5                                                            # E-value threshold for reporting hits

# Define desired output format columns for --outfmt 6
# Standard BLAST tabular fields plus query/subject lengths
# qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
OUTFMT_FIELDS="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# --- Functions ---
log_info() {
    echo "[INFO] $(date +'%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    # Print warning messages to stderr
    echo "[WARN] $(date +'%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_error() {
    # Print error messages to stderr
    echo "[ERROR] $(date +'%Y-%m-%d %H:%M:%S') - $1" >&2
}

# --- Pre-run Checks & Setup ---
log_info "Starting DIAMOND loop script."

# Check if diamond executable exists
if ! command -v diamond &> /dev/null; then
    log_error "DIAMOND executable not found in PATH. Please install DIAMOND (e.g., conda install diamond)."
    exit 1
fi
log_info "DIAMOND executable found: $(command -v diamond)"

# Check if target OG list exists
if [ ! -f "${TARGET_OG_LIST}" ]; then
    log_error "Target OG list file not found: ${TARGET_OG_LIST}"
    exit 1
fi
log_info "Using target OG list: ${TARGET_OG_LIST}"

# Check if input FASTA directory exists
if [ ! -d "${INPUT_FASTA_DIR}" ]; then
    log_error "Input FASTA directory not found: ${INPUT_FASTA_DIR}"
    exit 1
fi
log_info "Using input FASTA directory: ${INPUT_FASTA_DIR}"

# Check if DIAMOND database exists
DIAMOND_DB_FILE="${EUK_DB_NAME}.dmnd"
if [ ! -f "${DIAMOND_DB_FILE}" ]; then
    log_error "DIAMOND database file not found: ${DIAMOND_DB_FILE}"
    log_error "Ensure the database was created using 'diamond makedb'."
    # Example: diamond makedb --in ${EUK_DB_NAME}.fasta -d ${EUK_DB_NAME}
    exit 1
fi
log_info "Using DIAMOND database: ${DIAMOND_DB_FILE}"

# Create output directory
log_info "Ensuring output directory exists: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: ${OUTPUT_DIR}"
    exit 1
fi

# --- Main Loop ---
log_info "Starting DIAMOND searches..."
START_TIME=$SECONDS
COUNT=0
# Read target OGs into an array for potentially easier handling, though wc -l is fine too.
mapfile -t OGS_TO_PROCESS < <(grep -v '^\s*$' "${TARGET_OG_LIST}") # Read non-empty lines into array
TOTAL_OGS=${#OGS_TO_PROCESS[@]}
log_info "Found ${TOTAL_OGS} OGs to process."

PROCESSED_SUCCESSFULLY=0
SKIPPED_MISSING_FASTA=0
PROCESSED_WITH_ERRORS=0

# Iterate through the array of OG IDs
for OG_ID in "${OGS_TO_PROCESS[@]}"; do
    COUNT=$((COUNT + 1))
    QUERY_FASTA="${INPUT_FASTA_DIR}/${OG_ID}_phylo_input.fasta"
    OUTPUT_TSV="${OUTPUT_DIR}/${OG_ID}_diamond_hits.tsv"

    # Check if input FASTA exists for this OG
    if [ ! -f "${QUERY_FASTA}" ]; then
        log_warn "Input FASTA not found for OG ${COUNT}/${TOTAL_OGS} (${OG_ID}), skipping."
        SKIPPED_MISSING_FASTA=$((SKIPPED_MISSING_FASTA + 1))
        continue
    fi

    # Use printf for potentially cleaner, overwriting progress line if terminal supports it
    printf "[INFO] Processing OG %d/%d: %s ... \n" "$COUNT" "$TOTAL_OGS" "$OG_ID"

    # Run DIAMOND blastp
    # Use a temporary file for output, rename on success to avoid partial files on error
    TMP_OUTPUT_TSV="${OUTPUT_TSV}.tmp.$$"
    diamond blastp \
        --query "${QUERY_FASTA}" \
        --db "${DIAMOND_DB_FILE}" \
        --out "${TMP_OUTPUT_TSV}" \
        --outfmt 6 ${OUTFMT_FIELDS} \
        --sensitive \
        --max-target-seqs "${MAX_TARGETS}" \
        --evalue "${EVALUE}" \
        --threads "${THREADS}" \
        --quiet # Reduce verbosity per job

    # Check DIAMOND exit status
    DIAMOND_EXIT_CODE=$?
    if [ ${DIAMOND_EXIT_CODE} -ne 0 ]; then
        log_error "DIAMOND failed for ${OG_ID} (Exit code: ${DIAMOND_EXIT_CODE}). Check logs if available. Temp output: ${TMP_OUTPUT_TSV}"
        PROCESSED_WITH_ERRORS=$((PROCESSED_WITH_ERRORS + 1))
        # Optionally remove the temp file: rm -f "${TMP_OUTPUT_TSV}"
    else
        # Rename temp file to final name on success
        mv "${TMP_OUTPUT_TSV}" "${OUTPUT_TSV}"
        if [ $? -ne 0 ]; then
             log_error "Failed to rename temporary output file for ${OG_ID}."
             PROCESSED_WITH_ERRORS=$((PROCESSED_WITH_ERRORS + 1))
        else
            PROCESSED_SUCCESSFULLY=$((PROCESSED_SUCCESSFULLY + 1))
        fi
    fi
    # Clean up temp file if it still exists (e.g., if mv failed)
    rm -f "${TMP_OUTPUT_TSV}"

done

# --- Summary ---
ELAPSED_TIME=$((SECONDS - START_TIME))
log_info "--------------------------------------------------"
log_info "DIAMOND search loop finished."
log_info "Total OGs in list: ${TOTAL_OGS}"
log_info "  Processed successfully: ${PROCESSED_SUCCESSFULLY}"
log_info "  Skipped (missing FASTA): ${SKIPPED_MISSING_FASTA}"
log_info "  Failed (DIAMOND error): ${PROCESSED_WITH_ERRORS}"
log_info "Total execution time: ${ELAPSED_TIME} seconds."
log_info "Output TSV files are in: ${OUTPUT_DIR}"
log_info "--------------------------------------------------"

# Exit with error status if any DIAMOND runs failed
if [ ${PROCESSED_WITH_ERRORS} -gt 0 ]; then
    log_error "One or more DIAMOND runs failed. Please check logs."
    exit 1
else
    log_info "All DIAMOND runs completed successfully or were skipped as expected."
    exit 0
fi
