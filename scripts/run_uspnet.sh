#!/usr/bin/env bash

# Script to run USPNet-fast signal peptide prediction locally.
# v1.1 - Updated based on predict_fast.py script analysis.

# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# Prevent errors in a pipeline from being masked.
set -euo pipefail

# --- Default Configuration ---
DEFAULT_INPUT_FASTA="all_proteins_for_dtm.fasta"
# Directory to store intermediate processed data (relative to project root)
DEFAULT_PROCESSED_DATA_DIR="USPNet_Processed_Data"
# Directory to store final prediction output (relative to project root)
# Note: Final output CSV will be *inside* the processed data dir.
DEFAULT_OUTPUT_DIR_PARENT="USPNet_Results" # Parent dir for organization
# Path to the directory containing USPNet scripts (data_processing.py, predict_fast.py)
# Leave empty if scripts are directly in PATH after installation.
# If cloned, set to e.g., "path/to/cloned/USPNet"
DEFAULT_USPNET_DIR="" # MUST BE PROVIDED if not installed in PATH

# --- Usage Function ---
usage() {
  echo "Usage: $0 [OPTIONS]"
  echo "Runs USPNet-fast locally."
  echo "Ensure USPNet is installed (pip install or cloned) in the active Python environment."
  echo ""
  echo "Options:"
  echo "  -i, --input-fasta <file>    Path to the input FASTA file (Default: ${DEFAULT_INPUT_FASTA})"
  echo "  -p, --processed-dir <dir>   Directory for intermediate processed data (Default: ${DEFAULT_PROCESSED_DATA_DIR})"
  # echo "  -o, --output-dir <dir>      Parent directory for results (Unused, output goes to processed dir) (Default: ${DEFAULT_OUTPUT_DIR_PARENT})"
  echo "  -u, --uspnet-dir <dir>      Path to the USPNet script directory (REQUIRED if not installed in PATH) (Default: '${DEFAULT_USPNET_DIR}')"
  echo "  -h, --help                  Display this help message"
  echo ""
  echo "Example (if scripts in PATH - less likely for USPNet):"
  echo "  $0 -i all_proteins.fasta -p USPNet_Intermediate"
  echo "Example (if scripts in cloned dir):"
  echo "  $0 -i all_proteins.fasta -u path/to/USPNet -p USPNet_Intermediate"
  exit 1
}

# --- Parse Command-Line Arguments ---
INPUT_FASTA_ARG="$DEFAULT_INPUT_FASTA"
PROCESSED_DATA_DIR_ARG="$DEFAULT_PROCESSED_DATA_DIR"
# OUTPUT_DIR_PARENT_ARG="$DEFAULT_OUTPUT_DIR_PARENT" # Output dir arg removed
USPNET_DIR_ARG="$DEFAULT_USPNET_DIR"

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i|--input-fasta)
      INPUT_FASTA_ARG="$2"
      shift; shift
      ;;
    -p|--processed-dir)
      PROCESSED_DATA_DIR_ARG="$2"
      shift; shift
      ;;
    # -o|--output-dir) # Removed
    #   OUTPUT_DIR_PARENT_ARG="$2"
    #   shift; shift
    #   ;;
    -u|--uspnet-dir)
      USPNET_DIR_ARG="$2"
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

# Check if USPNet directory is provided if needed
if [ -z "$USPNET_DIR_ARG" ]; then
     if ! command -v data_processing.py &> /dev/null || ! command -v predict_fast.py &> /dev/null; then
          echo "ERROR: USPNet script directory (-u/--uspnet-dir) must be provided if scripts are not in PATH."
          usage
     fi
fi


echo "--- Starting USPNet-fast Prediction ---"

# Assumes running from the project root directory
PROJECT_DIR="$PWD"

# --- Construct Host Paths ---
INPUT_FASTA_HOST=$(cd "$(dirname "${PROJECT_DIR}/${INPUT_FASTA_ARG}")" && pwd)/$(basename "${PROJECT_DIR}/${INPUT_FASTA_ARG}")
PROCESSED_DATA_DIR_HOST="${PROJECT_DIR}/${PROCESSED_DATA_DIR_ARG}"
# OUTPUT_DIR_HOST="${PROJECT_DIR}/${OUTPUT_DIR_PARENT_ARG}" # Not directly used for output path

# Determine paths to USPNet scripts
if [ -n "$USPNET_DIR_ARG" ]; then
    # If USPNet dir provided, construct absolute path to scripts
    USPNET_DIR_ABS=$(cd "${PROJECT_DIR}/${USPNET_DIR_ARG}" && pwd) # Get absolute path to USPNet dir
    DATA_PROCESSING_SCRIPT="${USPNET_DIR_ABS}/data_processing.py"
    PREDICT_SCRIPT_NAME="predict_fast.py" # Just the name, will run from USPNet dir
    PREDICT_SCRIPT_PATH="${USPNET_DIR_ABS}/${PREDICT_SCRIPT_NAME}"
else
    # If USPNet dir not provided, assume scripts are in PATH
    DATA_PROCESSING_SCRIPT="data_processing.py"
    PREDICT_SCRIPT_NAME="predict_fast.py"
    PREDICT_SCRIPT_PATH=$(command -v ${PREDICT_SCRIPT_NAME}) # Find full path if in PATH
    USPNET_DIR_ABS=$(dirname "${PREDICT_SCRIPT_PATH}") # Infer dir from script path
fi


# --- Path/Command Validation ---
echo "--- Validating Paths & Commands ---"
if [ ! -f "$INPUT_FASTA_HOST" ]; then
    echo "ERROR: Input FASTA file not found: $INPUT_FASTA_HOST"
    exit 1
fi
echo "Input FASTA: $INPUT_FASTA_HOST"

# Check if scripts exist
if [ ! -f "$DATA_PROCESSING_SCRIPT" ]; then
    echo "ERROR: data_processing.py not found at expected location: $DATA_PROCESSING_SCRIPT"
    exit 1
fi
 if [ ! -f "$PREDICT_SCRIPT_PATH" ]; then
    echo "ERROR: predict_fast.py not found at expected location: $PREDICT_SCRIPT_PATH"
    exit 1
fi
echo "Using USPNet directory: $USPNET_DIR_ABS"
echo "Using data processing script: $DATA_PROCESSING_SCRIPT"
echo "Using prediction script: $PREDICT_SCRIPT_PATH"

# Check for model directory relative to USPNet dir
MODEL_DIR_REL="../data/mdl" # Relative path used in predict_fast.py
MODEL_DIR_ABS="${USPNET_DIR_ABS}/${MODEL_DIR_REL}"
if [ ! -d "$MODEL_DIR_ABS" ]; then
    echo "ERROR: Model directory not found at expected relative location: ${MODEL_DIR_ABS}"
    echo "       (Expected relative to USPNet script directory based on '../data/mdl' in predict_fast.py)"
    exit 1
fi
echo "Found model directory: ${MODEL_DIR_ABS}"


# Create directories
mkdir -p "$PROCESSED_DATA_DIR_HOST"
# mkdir -p "$OUTPUT_DIR_HOST" # Parent dir not strictly needed
echo "Intermediate data directory: $PROCESSED_DATA_DIR_HOST"
echo "Final output CSV will be inside: $PROCESSED_DATA_DIR_HOST"
echo "------------------------"

# --- Activate Conda Environment (Optional but Recommended) ---
# Ensure the environment containing USPNet and its dependencies is active.
# You can activate it manually before running this script, or uncomment below.
# echo "Activating conda environment 'tack_env_x86'..."
# source $(conda info --base)/etc/profile.d/conda.sh # Source conda functions
# conda activate tack_env_x86
# if [ $? -ne 0 ]; then
#     echo "ERROR: Failed to activate conda environment 'tack_env_x86'"
#     exit 1
# fi

# --- Step 1: Run Data Processing ---
echo
echo "Running USPNet Data Processing..."
DATA_PROC_CMD=(
    "python" "${DATA_PROCESSING_SCRIPT}"
      "--fasta_file" "${INPUT_FASTA_HOST}"
      "--data_processed_dir" "${PROCESSED_DATA_DIR_HOST}"
)
echo "--- Command ---"
printf "%q " "${DATA_PROC_CMD[@]}"
echo
echo "---------------"

if ! "${DATA_PROC_CMD[@]}"; then
    echo
    echo "--- USPNet Data Processing failed ---"
    # conda deactivate # Deactivate if activated by script
    exit 1
fi
echo "--- Data Processing finished successfully ---"


# --- Step 2: Run Prediction ---
echo
echo "Running USPNet Prediction (Fast)..."
# Run predict_fast.py from within its own directory to handle relative model paths
PREDICT_CMD=(
    "python" "${PREDICT_SCRIPT_NAME}" # Use just the script name relative to USPNet dir
      "--data_dir" "${PROCESSED_DATA_DIR_HOST}" # Use the processed data directory
      "--group_info" "no_group_info" # Use the no_group_info model
)
echo "--- Command (run from ${USPNET_DIR_ABS}) ---"
printf "%q " "${PREDICT_CMD[@]}"
echo
echo "------------------------------------------"

# Change directory, run prediction, change back
ORIGINAL_DIR="$PWD"
cd "$USPNET_DIR_ABS"

if ! "${PREDICT_CMD[@]}"; then
    EXIT_CODE=$?
    echo
    echo "--- USPNet Prediction failed with exit code $EXIT_CODE ---"
    cd "$ORIGINAL_DIR" # Change back even on failure
    # conda deactivate # Deactivate if activated by script
    exit $EXIT_CODE
fi

cd "$ORIGINAL_DIR" # Change back to original directory

echo
echo "--- USPNet Prediction completed successfully ---"
FINAL_OUTPUT_CSV="${PROCESSED_DATA_DIR_HOST}/results.csv"
echo "Output CSV file should be: ${FINAL_OUTPUT_CSV}"
if [ -f "$FINAL_OUTPUT_CSV" ]; then
    ls -l "$FINAL_OUTPUT_CSV"
else
    echo "WARNING: Expected output file ${FINAL_OUTPUT_CSV} not found!"
    echo "Listing contents of ${PROCESSED_DATA_DIR_HOST}:"
    ls -l "$PROCESSED_DATA_DIR_HOST"
fi


# Deactivate conda environment if activated by script
# conda deactivate

echo "Script finished."
exit 0
