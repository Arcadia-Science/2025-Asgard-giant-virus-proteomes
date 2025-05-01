#!/usr/bin/env python3

"""
Predicts intrinsic protein disorder using Metapredict and filters sequences
from input FASTA files based on the average disorder score and presence of
unknown ('X') amino acids.

Input: Directory containing FASTA files (e.g., length-filtered).
Output:
- Directory containing FASTA files of predicted 'globular' proteins (below threshold, no 'X').
- FASTA file containing sequences skipped due to 'X'.
- FASTA file containing sequences predicted as 'disordered' (at or above threshold, no 'X').
"""

import os
import glob
import sys
import numpy as np
from Bio import SeqIO
import argparse # For command-line arguments
import time # For timing

# --- Library Imports ---
try:
    import metapredict
except ImportError:
    print("ERROR: Metapredict library not found. Please install it (`pip install metapredict`).", file=sys.stderr)
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    print("ERROR: Biopython library not found. Please install it (`pip install biopython`).", file=sys.stderr)
    sys.exit(1)

# --- Constants ---
# Default values - can be overridden by command-line arguments
DEFAULT_INPUT_DIR = 'Asgard_fasta_files_length_filtered'
DEFAULT_GLOBULAR_DIR = 'Asgard_faa_globular_filtered'
DEFAULT_SKIPPED_X_FILE = 'Asgard_proteins_with_X.fasta'
DEFAULT_DISORDERED_FILE = 'Asgard_proteins_predicted_disordered.fasta'
DEFAULT_DISORDER_THRESHOLD = 0.5 # Score >= threshold -> Disordered

# --- Function Definitions ---

def setup_logging(level='INFO'):
    """Basic logging setup."""
    # In a larger script, consider using the 'logging' module
    def log(msg_level, message):
        if level == 'INFO' or msg_level == 'ERROR' or msg_level == 'WARN':
             # Simple timestamp
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            output_stream = sys.stderr if msg_level in ['ERROR', 'WARN'] else sys.stdout
            print(f"[{msg_level}] {timestamp} - {message}", file=output_stream)
    return log

def setup_output_dir(dir_path, log):
    """Creates an output directory if it doesn't exist."""
    try:
        os.makedirs(dir_path, exist_ok=True)
        log('INFO', f"Output directory ensured: {dir_path}")
    except OSError as e:
        log('ERROR', f"Creating output directory {dir_path}: {e}")
        sys.exit(1)

def process_fasta_file(input_filepath, globular_dir, skipped_writer, disordered_writer, threshold, log):
    """
    Processes a single FASTA file, predicts disorder, filters, and writes outputs
    using the provided writer functions/handles.

    Args:
        input_filepath (str): Path to the input FASTA file.
        globular_dir (str): Directory to write globular output FASTA.
        skipped_writer (function/None): Function to write skipped ('X') records or None.
        disordered_writer (function/None): Function to write disordered records or None.
        threshold (float): Disorder threshold.
        log (function): Logging function.

    Returns:
        tuple: (records_in_file, records_kept_globular, records_skipped_x, records_discarded_disordered)
    """
    basename = os.path.basename(input_filepath)
    globular_output_filepath = os.path.join(globular_dir, basename)

    # Counters for this specific file
    file_total = 0
    file_globular = 0
    file_skipped_x = 0
    file_disordered = 0
    file_error_predict = 0
    file_error_other = 0

    globular_records_to_write = [] # Collect records to write at the end for this file

    try:
        for record in SeqIO.parse(input_filepath, "fasta"):
            file_total += 1

            # --- Initial Sequence Validation ---
            if not record.seq:
                log('WARN', f"Skipping record '{record.id}' in '{basename}' due to empty sequence.")
                file_error_other += 1
                continue
            sequence = str(record.seq).upper()
            if not sequence: # Check again after potential conversion issues
                log('WARN', f"Skipping record '{record.id}' in '{basename}' (empty after conversion).")
                file_error_other += 1
                continue

            # --- Check for 'X' ---
            if 'X' in sequence:
                file_skipped_x += 1
                if skipped_writer:
                    skipped_writer(record)
                continue # Skip prediction for sequences with 'X'

            # --- Predict Disorder ---
            try:
                # Using CPU explicitly as per original script's intent
                disorder_scores = metapredict.predict_disorder(sequence, device='cpu')
            except Exception as predict_err:
                log('ERROR', f"Metapredict failed for record '{record.id}' in '{basename}': {predict_err}")
                file_error_predict += 1
                continue # Skip this record

            # --- Process Scores ---
            if disorder_scores is None:
                 log('WARN', f"Metapredict returned None for record '{record.id}' in '{basename}'. Skipping.")
                 file_error_predict += 1
                 continue
            try:
                average_disorder = np.mean(disorder_scores)
                # Check for NaN (Not a Number) result
                if np.isnan(average_disorder):
                    log('WARN', f"Metapredict average score is NaN for record '{record.id}' in '{basename}'. Skipping.")
                    file_error_predict += 1
                    continue
            except Exception as mean_err:
                 log('ERROR', f"Calculating mean disorder for '{record.id}' in '{basename}': {mean_err}")
                 file_error_other += 1
                 continue

            # --- Apply Filter ---
            if average_disorder < threshold:
                # Keep as globular
                globular_records_to_write.append(record)
                file_globular += 1
            else:
                # Discard as disordered
                file_disordered += 1
                if disordered_writer:
                    disordered_writer(record)

        # --- Write Globular Output for this File ---
        if globular_records_to_write:
            try:
                # Writing all at once is generally efficient
                with open(globular_output_filepath, "w") as outfile:
                    SeqIO.write(globular_records_to_write, outfile, "fasta")
            except IOError as write_err:
                 log('ERROR', f"Failed to write globular output file '{globular_output_filepath}': {write_err}")
                 # Note: Some records for this file might be lost if writing fails here.

    except Exception as e:
        # Catch unexpected errors during file parsing or processing
        log('ERROR', f"Unhandled exception processing file '{basename}': {e}")
        # Re-raise to stop execution, as state might be inconsistent
        raise

    # Return counts for this file
    return file_total, file_globular, file_skipped_x, file_disordered, file_error_predict, file_error_other


def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Predict protein disorder and filter FASTA files.")
    parser.add_argument(
        "-i", "--input-dir",
        default=DEFAULT_INPUT_DIR,
        help=f"Input directory containing FASTA files (default: {DEFAULT_INPUT_DIR})"
    )
    parser.add_argument(
        "-g", "--globular-dir",
        default=DEFAULT_GLOBULAR_DIR,
        help=f"Output directory for globular proteins (default: {DEFAULT_GLOBULAR_DIR})"
    )
    parser.add_argument(
        "-x", "--skipped-file",
        default=DEFAULT_SKIPPED_X_FILE,
        help=f"Output FASTA file for sequences skipped due to 'X' (default: {DEFAULT_SKIPPED_X_FILE})"
    )
    parser.add_argument(
        "-d", "--disordered-file",
        default=DEFAULT_DISORDERED_FILE,
        help=f"Output FASTA file for predicted disordered sequences (default: {DEFAULT_DISORDERED_FILE})"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=DEFAULT_DISORDER_THRESHOLD,
        help=f"Disorder threshold (score >= threshold is disordered) (default: {DEFAULT_DISORDER_THRESHOLD})"
    )
    return parser.parse_args()

def main():
    """Main execution logic."""
    # --- Argument Parsing & Setup ---
    args = parse_arguments()
    log = setup_logging() # Initialize logging

    log('INFO', "Starting Disorder Prediction and Filtering script.")
    log('INFO', f"Input FASTA directory: {args.input_dir}")
    log('INFO', f"Output globular directory: {args.globular_dir}")
    log('INFO', f"Output skipped ('X') file: {args.skipped_file}")
    log('INFO', f"Output disordered file: {args.disordered_file}")
    log('INFO', f"Disorder threshold: >= {args.threshold}")

    # Validate Input Directory
    if not os.path.isdir(args.input_dir):
        log('ERROR', f"Input directory '{args.input_dir}' not found.")
        sys.exit(1)

    # Setup Output Directory
    setup_output_dir(args.globular_dir, log)

    # Find Input FASTA Files
    input_fasta_files = glob.glob(os.path.join(args.input_dir, '*.fasta')) + \
                        glob.glob(os.path.join(args.input_dir, '*.faa')) # Add .faa common extension
    if not input_fasta_files:
        log('ERROR', f"No '.fasta' or '.faa' files found in '{args.input_dir}'.")
        sys.exit(1)
    log('INFO', f"Found {len(input_fasta_files)} input FASTA files to process.")

    # --- Initialize Counters ---
    total_records_processed = 0
    total_records_kept_globular = 0
    total_files_processed = 0
    total_skipped_x = 0
    total_discarded_disordered = 0
    total_prediction_errors = 0
    total_other_errors = 0
    failed_files = []

    # --- Process Files ---
    log('INFO', "Starting prediction and filtering loop...")
    start_time = time.time()

    # Use 'with open' for the summary files to ensure they are closed
    try:
        with open(args.skipped_file, "w") as skipped_fh, \
             open(args.disordered_file, "w") as disordered_fh:

            # Define writer functions to pass to the processing function
            def skipped_writer(record):
                SeqIO.write([record], skipped_fh, "fasta")

            def disordered_writer(record):
                SeqIO.write([record], disordered_fh, "fasta")

            for i, input_filepath in enumerate(sorted(input_fasta_files)):
                basename = os.path.basename(input_filepath)
                # Use carriage return for progress update on the same line
                print(f"  Processing file {i+1}/{len(input_fasta_files)}: {basename}...", end='\r', flush=True)

                try:
                    file_processed, file_kept, file_skipped, file_discarded, file_err_pred, file_err_other = process_fasta_file(
                        input_filepath,
                        args.globular_dir,
                        skipped_writer,
                        disordered_writer,
                        args.threshold,
                        log # Pass the log function
                    )
                    # Accumulate counts
                    total_files_processed += 1
                    total_records_processed += file_processed
                    total_records_kept_globular += file_kept
                    total_skipped_x += file_skipped
                    total_discarded_disordered += file_discarded
                    total_prediction_errors += file_err_pred
                    total_other_errors += file_err_other

                except Exception as proc_err:
                    # Catch errors raised by process_fasta_file (should be rare now)
                    log('ERROR', f"FATAL unhandled exception during processing of {basename}: {proc_err}")
                    failed_files.append(basename)
                    # Decide whether to stop or continue
                    # For now, log and continue with the next file
                    log('WARN', "Attempting to continue with the next file...")

        # Clear the progress line after the loop
        print(" " * 80, end='\r')

    except IOError as e:
        log('ERROR', f"Could not open summary output file {e.filename} for writing: {e}")
        # Can't proceed without the summary files in this setup
        sys.exit(1)

    end_time = time.time()
    elapsed_time = end_time - start_time

    # --- Final Summary ---
    log('INFO', "--------------------------------------------------")
    log('INFO', "Disorder prediction and filtering complete.")
    log('INFO', f"Processed {total_files_processed}/{len(input_fasta_files)} input files in {elapsed_time:.2f} seconds.")
    if failed_files:
        log('WARN', f"Processing failed for {len(failed_files)} files: {', '.join(failed_files)}")
    log('INFO', f"Total proteins encountered: {total_records_processed:,}")
    log('INFO', f"Total proteins skipped due to 'X': {total_skipped_x:,}")
    if total_skipped_x > 0:
        log('INFO', f"Skipped ('X') sequences saved to: {args.skipped_file}")
    log('INFO', f"Total predicted disordered (score >= {args.threshold}, no 'X'): {total_discarded_disordered:,}")
    if total_discarded_disordered > 0:
         log('INFO', f"Predicted disordered sequences saved to: {args.disordered_file}")
    log('INFO', f"Total predicted globular (score < {args.threshold}, no 'X'): {total_records_kept_globular:,}")
    log('INFO', f"Total Metapredict errors encountered: {total_prediction_errors:,}")
    log('INFO', f"Total other processing errors (e.g., empty seq): {total_other_errors:,}")


    # Calculate percentage based on proteins actually submitted for prediction
    valid_for_prediction = total_records_processed - total_skipped_x - total_prediction_errors - total_other_errors
    if valid_for_prediction > 0:
        kept_percentage = (total_records_kept_globular / valid_for_prediction) * 100
        discarded_percentage = (total_discarded_disordered / valid_for_prediction) * 100
        log('INFO', f"Of {valid_for_prediction:,} proteins submitted for prediction (no 'X', no errors):")
        log('INFO', f"  - Kept as globular: {total_records_kept_globular:,} ({kept_percentage:.1f}%)")
        log('INFO', f"  - Discarded as disordered: {total_discarded_disordered:,} ({discarded_percentage:.1f}%)")
    elif total_records_processed > 0:
         log('WARN', "No valid proteins were successfully submitted for prediction (check for 'X' or errors).")
    else:
        log('INFO', "No proteins found in input files.")

    log('INFO', f"Filtered globular FASTA files are in: {args.globular_dir}")
    log('INFO', "--------------------------------------------------")

    # Reminder for next steps
    print(f"\nNext step suggestion:")
    print(f"Run OrthoFinder using '{args.globular_dir}' as the input directory ('-f' option).")

    if total_prediction_errors > 0 or total_other_errors > 0 or failed_files:
        log('WARN', "Script finished with some errors or warnings. Please review the log.")
        sys.exit(1) # Exit with non-zero status if errors occurred
    else:
        log('INFO', "Script finished successfully.")
        sys.exit(0)


if __name__ == "__main__":
    main()

