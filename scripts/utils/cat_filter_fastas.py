#!/usr/bin/env python3

"""
Prepares FASTA subsets for downstream analysis (e.g., InterProScan, CD-HIT).

1. Concatenates multiple FASTA files from an input directory into a single file.
2. Filters the concatenated FASTA file based on keywords in the sequence headers
   to create a subset (e.g., hypothetical proteins).
"""

import os
import glob
import sys
import re
import argparse
import logging
from pathlib import Path
from Bio import SeqIO # Requires biopython

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# --- Default Keywords for Filtering (if extracting hypotheticals) ---
# These can be overridden via command-line arguments.
DEFAULT_KEYWORDS_TO_KEEP = [
    "hypothetical",
    "uncharacterized",
    "unknown function",
    "predicted protein",
    "conserved protein",
    "putative protein",
    "DUF", # Domain of Unknown Function
    "unnamed protein product",
    "orf", # Open reading frame - use with caution
    "possible protein"
]

# --- Function Definitions ---

def concatenate_fastas(input_dir, output_file):
    """
    Concatenates all FASTA files found in a directory into a single output file.

    Args:
        input_dir (str): Path to the directory containing input FASTA files.
                         Files ending with '.fasta', '.faa', '.fa' are considered.
        output_file (str): Path to the single output FASTA file to be created.

    Returns:
        bool: True if successful, False otherwise.
    """
    logging.info(f"Starting FASTA concatenation from '{input_dir}' to '{output_file}'...")
    input_path = Path(input_dir)
    output_path = Path(output_file)
    fasta_files = []

    if not input_path.is_dir():
        logging.error(f"Input directory not found: {input_dir}")
        return False

    # Ensure the output directory exists
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logging.error(f"Cannot create output directory {output_path.parent}: {e}")
        return False

    # Find all fasta files (common extensions)
    for ext in ['*.fasta', '*.faa', '*.fa']:
        fasta_files.extend(input_path.glob(ext))

    # Remove duplicates and sort
    fasta_files = sorted(list(set(fasta_files)))

    if not fasta_files:
        logging.error(f"No FASTA files (.fasta, .faa, .fa) found in directory: {input_dir}")
        return False

    logging.info(f"Found {len(fasta_files)} FASTA files to concatenate.")

    try:
        total_records_written = 0
        with open(output_path, 'w') as outfile:
            processed_count = 0
            for i, fname in enumerate(fasta_files):
                # Simple progress indicator
                if (i + 1) % 50 == 0 or (i + 1) == len(fasta_files):
                     logging.info(f"  Processing file {i+1}/{len(fasta_files)}: {fname.name}")
                try:
                    records_in_file = 0
                    with open(fname, 'r') as infile:
                        for record in SeqIO.parse(infile, "fasta"):
                            # Basic validation
                            if not record.id or not record.seq:
                                logging.warning(f"  Skipping potentially invalid record in {fname.name} (ID: {record.id[:50]}...)")
                                continue
                            SeqIO.write(record, outfile, "fasta")
                            records_in_file += 1
                    total_records_written += records_in_file
                    processed_count += 1
                except Exception as e:
                    logging.warning(f"  Could not fully read or process file {fname.name}. Skipping remaining. Error: {e}")
                    continue # Skip to next file

        logging.info(f"Concatenation complete.")
        logging.info(f"Processed {processed_count}/{len(fasta_files)} files.")
        logging.info(f"Wrote {total_records_written:,} total records to: {output_file}")
        return True

    except IOError as e:
        logging.error(f"Error writing to output file {output_file}: {e}")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred during concatenation: {e}")
        return False


def filter_fasta_by_keywords(input_fasta_path, output_fasta_path, keywords_to_keep, header_delimiter='|', name_field_index=4):
    """
    Reads an input FASTA file and writes sequences to an output FASTA file
    if their header matches specific criteria based on keywords in a designated field.

    Args:
        input_fasta_path (str): Path to the input concatenated FASTA file.
        output_fasta_path (str): Path to write the filtered output FASTA file.
        keywords_to_keep (list): A list of keywords (strings) to search for
                                 (case-insensitive) in the designated name field.
                                 Sequences are kept if the name field is blank OR
                                 contains any of these keywords.
        header_delimiter (str): The delimiter used in the FASTA header. Default '|'.
        name_field_index (int): The 0-based index of the field containing the protein
                                name/description after splitting by the delimiter. Default 4.

    Returns:
        bool: True if successful, False otherwise.
    """
    logging.info(f"Starting FASTA filtering from '{input_fasta_path}' to '{output_fasta_path}'...")
    logging.info(f"  Keeping sequences if field {name_field_index+1} (delimited by '{header_delimiter}') is blank OR contains (case-insensitive): {', '.join(keywords_to_keep)}")

    input_path = Path(input_fasta_path)
    output_path = Path(output_fasta_path)

    if not input_path.is_file():
        logging.error(f"Input file not found: {input_path}")
        return False

    # Ensure output directory exists
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logging.error(f"Cannot create output directory {output_path.parent}: {e}")
        return False

    # Compile regex patterns for efficiency (case-insensitive)
    patterns = [re.compile(re.escape(kw), re.IGNORECASE) for kw in keywords_to_keep]

    sequences_written = 0
    sequences_read = 0
    try:
        infile_size = input_path.stat().st_size
        logging.info(f"Input file size: {infile_size / (1024*1024):.2f} MB")

        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                sequences_read += 1
                if sequences_read % 50000 == 0:
                    logging.info(f"  Processed {sequences_read:,} headers...")

                keep_record = False
                name_part = ""
                try:
                    # Use record.id for the part before the first space
                    header_id = record.id
                    parts = header_id.split(header_delimiter)

                    if len(parts) > name_field_index:
                        name_part = parts[name_field_index].strip()
                        if not name_part: # Keep if name field is blank
                            keep_record = True
                        else:
                            # Check if any keyword matches
                            for pattern in patterns:
                                if pattern.search(name_part):
                                    keep_record = True
                                    break # Found a match, no need to check further keywords
                    else:
                        # Header format is unexpected, treat as potentially unknown/hypothetical
                        logging.warning(f"  Header format incorrect (less than {name_field_index+1} parts delimited by '{header_delimiter}'): {header_id[:100]}...")
                        keep_record = True # Keep records with unexpected headers by default

                except Exception as e:
                    logging.warning(f"  Error parsing header '{record.id[:100]}...': {e}. Treating as keep.")
                    keep_record = True # Keep record if header parsing fails

                if keep_record:
                    SeqIO.write(record, outfile, "fasta")
                    sequences_written += 1

        logging.info(f"Filtering complete.")
        logging.info(f"Read {sequences_read:,} total sequences.")
        logging.info(f"Wrote {sequences_written:,} sequences matching criteria to: {output_path}")
        return True

    except FileNotFoundError:
        logging.error(f"Input file not found during processing: {input_path}")
        return False
    except IOError as e:
        logging.error(f"Error writing to output file {output_path}: {e}")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred during filtering: {e}")
        return False


def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Concatenate and filter FASTA files based on header keywords.")

    parser.add_argument("-i", "--input_dir", required=True,
                        help="Input directory containing individual FASTA files (.fa, .faa, .fasta).")
    parser.add_argument("-c", "--concat_output", required=True,
                        help="Output file path for the concatenated FASTA.")
    parser.add_argument("-f", "--filter_output", required=True,
                        help="Output file path for the filtered FASTA subset.")
    parser.add_argument("-k", "--keywords", nargs='+', default=DEFAULT_KEYWORDS_TO_KEEP,
                        help="List of keywords (case-insensitive) to keep in the filtered subset. "
                             "Sequences are kept if the name field is blank OR contains any keyword. "
                             f"(Default: {' '.join(DEFAULT_KEYWORDS_TO_KEEP)})")
    parser.add_argument("--delimiter", default='|',
                        help="Delimiter used in FASTA headers to separate fields (default: '|').")
    parser.add_argument("--name_index", type=int, default=4,
                        help="0-based index of the field containing the name/description "
                             "after splitting the header by the delimiter (default: 4).")

    return parser.parse_args()

# --- Main Execution ---
if __name__ == "__main__":
    args = parse_arguments()

    # Step 1: Concatenate
    if concatenate_fastas(args.input_dir, args.concat_output):
        # Step 2: Filter (only if concatenation succeeded)
        if not filter_fasta_by_keywords(args.concat_output, args.filter_output, args.keywords, args.delimiter, args.name_index):
            logging.error("FASTA filtering step failed.")
            sys.exit(1)
    else:
        logging.error("FASTA concatenation step failed.")
        sys.exit(1)

    logging.info("Script finished successfully.")
