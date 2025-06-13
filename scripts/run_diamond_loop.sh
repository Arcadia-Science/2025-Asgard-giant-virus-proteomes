#!/bin/bash

# ==============================================================================
# This script runs DIAMOND blastp for all .faa files in a given input directory
# and saves the results to a specified output directory.
#
# Usage: ./run_diamond_loop.sh <input_directory> <output_directory>
#   - input_directory:  Path to the directory containing .faa protein files.
#   - output_directory: Path to the directory where DIAMOND results will be saved.
# ==============================================================================

# --- ADDED: Argument parsing and validation ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
DB_PATH="../data/reference/uniref100.dmnd" # Assumes a fixed location for the database

# --- ADDED: Create output directory if it doesn't exist ---
mkdir -p "$OUTPUT_DIR"

# Check if the database file exists
if [ ! -f "$DB_PATH" ]; then
    echo "Error: DIAMOND database not found at $DB_PATH"
    exit 1
fi

echo "Starting DIAMOND blastp searches..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Loop through all .faa files in the input directory
for f in "$INPUT_DIR"/*.faa; do
    # Get the base name of the file
    bname=$(basename "$f" .faa)
    echo "Processing $bname"

    # Define the output file path
    output_file="$OUTPUT_DIR/${bname}_diamond_hits.tsv"

    # Run diamond blastp
    diamond blastp \
        -d "$DB_PATH" \
        -q "$f" \
        -o "$output_file" \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --threads 16 \
        --max-target-seqs 5 \
        --sensitive
done

echo "DIAMOND runs completed."