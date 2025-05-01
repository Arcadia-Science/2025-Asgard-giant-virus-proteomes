#!/usr/bin/env python3

"""
Filters sequences in FASTA files based on minimum and maximum length thresholds.

Reads FASTA files from an input directory, filters each sequence based on the
provided length criteria, and writes the kept sequences to corresponding
new FASTA files in an output directory.
"""

import os
import sys
import argparse
import logging
from pathlib import Path
import glob
from Bio import SeqIO # Requires biopython

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def filter_single_fasta(input_filepath: Path, output_filepath: Path, min_len: int, max_len: int) -> tuple[int, int]:
    """
    Filters a single FASTA file by sequence length.

    Args:
        input_filepath: Path to the input FASTA file.
        output_filepath: Path to write the filtered output FASTA file.
        min_len: Minimum sequence length (inclusive).
        max_len: Maximum sequence length (inclusive).

    Returns:
        Tuple (records_read, records_written).
    """
    records_read = 0
    records_written = 0
    filtered_records = []

    try:
        for record in SeqIO.parse(input_filepath, "fasta"):
            records_read += 1
            seq_len = len(record.seq)

            # Apply length filter
            if min_len <= seq_len <= max_len:
                filtered_records.append(record)
                records_written += 1

        # Write the filtered records to the new file, ONLY if any were kept
        if filtered_records:
            # Ensure output directory exists just before writing
            output_filepath.parent.mkdir(parents=True, exist_ok=True)
            with open(output_filepath, 'w') as outfile:
                SeqIO.write(filtered_records, outfile, "fasta")
        else:
            logging.debug(f"  No records kept for {input_filepath.name} after length filtering.")

        return records_read, records_written

    except FileNotFoundError:
        logging.error(f"Input file not found during processing: {input_filepath}")
        return 0, 0 # Indicate file couldn't be processed
    except Exception as e:
        logging.error(f"Error processing file {input_filepath.name}: {e}")
        return records_read, 0 # Return read count, but 0 written due to error

def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Filter FASTA files by sequence length.")
    parser.add_argument("-i", "--input_dir", required=True, type=Path,
                        help="Input directory containing FASTA files (.fa, .faa, .fasta).")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Output directory for the length-filtered FASTA files.")
    parser.add_argument("--min_len", type=int, default=50,
                        help="Minimum sequence length to keep (inclusive, default: 50).")
    parser.add_argument("--max_len", type=int, default=100000, # Set a high default max
                        help="Maximum sequence length to keep (inclusive, default: 100000).")
    return parser.parse_args()

# --- Main Execution ---
if __name__ == "__main__":
    args = parse_arguments()

    logging.info("Starting FASTA length filtering script...")
    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Length filter: {args.min_len} <= length <= {args.max_len}")

    # --- Validate Input / Setup Output ---
    if not args.input_dir.is_dir():
        logging.error(f"Input directory '{args.input_dir}' not found.")
        sys.exit(1)
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Output directory ensured: {args.output_dir}")
    except OSError as e:
        logging.error(f"Error creating output directory {args.output_dir}: {e}")
        sys.exit(1)

    # --- Find Input FASTA Files ---
    input_fasta_files = list(args.input_dir.glob('*.fasta')) + \
                        list(args.input_dir.glob('*.faa')) + \
                        list(args.input_dir.glob('*.fa'))
    input_fasta_files = sorted(list(set(input_fasta_files))) # Unique and sorted

    if not input_fasta_files:
        logging.error(f"No FASTA files (.fasta, .faa, .fa) found in '{args.input_dir}'.")
        sys.exit(1)
    logging.info(f"Found {len(input_fasta_files)} FASTA files to process.")

    # --- Process and Filter Each File ---
    total_records_processed = 0
    total_records_kept = 0
    total_files_processed = 0
    files_with_errors = 0

    for i, input_filepath in enumerate(input_fasta_files):
        total_files_processed += 1
        output_filepath = args.output_dir / input_filepath.name # Keep original filename
        # Use carriage return for cleaner progress updates
        print(f"  Processing file {total_files_processed}/{len(input_fasta_files)}: {input_filepath.name}...", end='\\r', flush=True)

        read_count, written_count = filter_single_fasta(input_filepath, output_filepath, args.min_len, args.max_len)

        total_records_processed += read_count
        total_records_kept += written_count
        if read_count > 0 and written_count == 0 and not Path(output_filepath).exists():
             # Check if an error occurred vs. just no sequences passing filter
             # This check is imperfect but better than nothing
             files_with_errors +=1


    # --- Final Summary ---
    print(" " * 80, end='\\r') # Clear the progress line
    logging.info("="*60)
    logging.info("Length filtering complete.")
    logging.info(f"Processed {total_files_processed} input files.")
    logging.info(f"Total sequences read: {total_records_processed:,}")
    logging.info(f"Total sequences kept (length {args.min_len}-{args.max_len} AA): {total_records_kept:,}")
    if total_records_processed > 0:
        kept_percentage = (total_records_kept / total_records_processed) * 100
        logging.info(f"Kept approximately {kept_percentage:.1f}% of the sequences.")
    else:
        logging.info("No sequences processed.")
    if files_with_errors > 0:
         logging.warning(f"Encountered errors while processing {files_with_errors} files. Check logs above.")
    logging.info(f"Filtered FASTA files are in: {args.output_dir}")
    logging.info("="*60)

    if files_with_errors > 0:
         sys.exit(1) # Exit with error if problems occurred
