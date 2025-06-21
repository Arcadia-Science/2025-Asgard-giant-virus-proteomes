#!/bin/bash

# ==============================================================================
# This script runs DIAMOND blastp for all protein FASTA files (e.g., .faa, .fasta)
# in a given input directory and saves the results to a specified output directory.
#
# It allows configuration of the DIAMOND database, number of threads, and
# max target sequences.
# ==============================================================================

set -euo pipefail

# --- Default Configuration ---
DEFAULT_DB_PATH="../data/reference/uniref100.dmnd"
DEFAULT_THREADS="16"
DEFAULT_MAX_TARGET_SEQS="5"
DEFAULT_INPUT_PATTERN="*.faa" # Default to .faa, can be changed

# --- Usage Function ---
usage() {
  echo "Usage: $0 -i <input_directory> -o <output_directory> [OPTIONS]"
  echo ""
  echo "Required arguments:"
  echo "  -i, --input-dir <dir>     Path to the directory containing protein FASTA files."
  echo "  -o, --output-dir <dir>    Path to the directory where DIAMOND results will be saved."
  echo ""
  echo "Optional arguments:"
  echo "  -d, --db <path>           Path to the DIAMOND database (default: ${DEFAULT_DB_PATH})"
  echo "  -t, --threads <int>       Number of threads for DIAMOND (default: ${DEFAULT_THREADS})"
  echo "  -m, --max-seqs <int>      Maximum number of target sequences per query for DIAMOND (default: ${DEFAULT_MAX_TARGET_SEQS})"
  echo "  -p, --pattern <glob>      Glob pattern for input FASTA files (default: '${DEFAULT_INPUT_PATTERN}') e.g., '*.fasta'"
  echo "  -h, --help                Display this help message."
  exit 1
}

# --- Parse Command-Line Arguments ---
INPUT_DIR_ARG=""
OUTPUT_DIR_ARG=""
DB_PATH_ARG="$DEFAULT_DB_PATH"
THREADS_ARG="$DEFAULT_THREADS"
MAX_TARGET_SEQS_ARG="$DEFAULT_MAX_TARGET_SEQS"
INPUT_PATTERN_ARG="$DEFAULT_INPUT_PATTERN"

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i|--input-dir)
      INPUT_DIR_ARG="$2"
      shift; shift
      ;;
    -o|--output-dir)
      OUTPUT_DIR_ARG="$2"
      shift; shift
      ;;
    -d|--db)
      DB_PATH_ARG="$2"
      shift; shift
      ;;
    -t|--threads)
      THREADS_ARG="$2"
      shift; shift
      ;;
    -m|--max-seqs)
      MAX_TARGET_SEQS_ARG="$2"
      shift; shift
      ;;
    -p|--pattern)
      INPUT_PATTERN_ARG="$2"
      shift; shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "ERROR: Unknown option '$1'"
      usage
      ;;
  esac
done

# Validate required arguments
if [ -z "$INPUT_DIR_ARG" ] || [ -z "$OUTPUT_DIR_ARG" ]; then
    echo "ERROR: Input and output directories (-i, -o) are required."
    usage
fi

# --- Setup ---
# Use absolute paths for robustness
INPUT_DIR="$(cd "$INPUT_DIR_ARG" && pwd)"
OUTPUT_DIR="$(cd "$OUTPUT_DIR_ARG" && pwd)" # Will be created if it doesn't exist

# If DB_PATH_ARG is relative, make it relative to PWD from where script is called, or assume absolute
# For simplicity, this example assumes user provides a usable path (absolute or correctly relative)
DB_PATH="$DB_PATH_ARG"

mkdir -p "$OUTPUT_DIR"

if [ ! -f "$DB_PATH" ]; then
    echo "ERROR: DIAMOND database not found at '$DB_PATH'"
    exit 1
fi

if ! command -v diamond &> /dev/null; then
    echo "ERROR: diamond executable not found in PATH. Please install DIAMOND."
    exit 1
fi

echo "--- Starting DIAMOND blastp searches ---"
echo "Input directory:   $INPUT_DIR"
echo "Input pattern:     $INPUT_PATTERN_ARG"
echo "Output directory:  $OUTPUT_DIR"
echo "DIAMOND DB:        $DB_PATH"
echo "Threads:           $THREADS_ARG"
echo "Max target seqs:   $MAX_TARGET_SEQS_ARG"
echo "----------------------------------------"

shopt -s nullglob # Important: ensures loop doesn't run if no files match

file_count=0
processed_count=0

for f in "$INPUT_DIR"/$INPUT_PATTERN_ARG; do
    file_count=$((file_count + 1))
    # Get the base name of the file, removing common FASTA extensions
    # This is a bit more robust for various extensions like .fasta, .faa, .fa
    bname=$(basename "$f")
    bname_noext="${bname%.*}" # Removes last extension
    if [[ "$bname_noext" == "$bname" ]]; then # If no extension, use full name (should not happen with pattern)
      bname_noext=$bname
    fi
    # Handle potential double extensions like .fasta.gz by removing common ones iteratively
    # For this script, assuming simple extensions as per pattern.

    echo "Processing $bname_noext ($f)..."

    output_file="$OUTPUT_DIR/${bname_noext}_diamond_hits.tsv"

    # Run diamond blastp
    if diamond blastp \
        --db "$DB_PATH" \
        --query "$f" \
        --out "$output_file" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --threads "$THREADS_ARG" \
        --max-target-seqs "$MAX_TARGET_SEQS_ARG" \
        --sensitive; then
      echo "Successfully processed $bname_noext."
      processed_count=$((processed_count + 1))
    else
      echo "ERROR: DIAMOND failed for $bname_noext. Check output/errors above."
      # Continue to next file as per original behavior, set -e handles script exit on diamond error
    fi
done

shopt -u nullglob # Reset nullglob

echo "--- DIAMOND runs completed ---"
if [ "$file_count" -eq 0 ]; then
  echo "No files found matching pattern '$INPUT_PATTERN_ARG' in '$INPUT_DIR'."
else
  echo "Processed $processed_count out of $file_count files found."
fi

if [ "$processed_count" -lt "$file_count" ] && [ "$file_count" -gt 0 ]; then
    echo "WARNING: Some files may have failed processing."
    exit 1 # Exit with error if not all files processed successfully
fi

exit 0