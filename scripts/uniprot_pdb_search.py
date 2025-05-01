#!/usr/bin/env python3

"""
Queries the UniProt REST API for individual UniProtKB accessions to retrieve
associated PDB cross-references.

Reads a list of UniProtKB ACs from an input file, queries the API for each,
and saves the results (UniProtKB_AC, PDB_IDs) to a TSV file.
Supports resuming from a previous run by checking the output file.
"""

import requests
import time
import os
import io # To handle string IO for TSV parsing
import argparse
import logging
import csv
from pathlib import Path
import sys
import shutil # To check for requests library

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# --- Constants ---
API_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}"
FIELDS_TO_RETRIEVE = "xref_pdb" # Field containing PDB cross-references
RESPONSE_FORMAT = "tsv"

# --- Helper Functions ---

def read_uniprot_ids(filename: Path) -> list[str]:
    """Reads unique UniProt IDs from a file, one ID per line."""
    if not filename.exists():
        logging.error(f"Input file not found at {filename}")
        raise FileNotFoundError(f"Input file not found: {filename}")
    try:
        with open(filename, 'r') as f:
            # Read, strip whitespace, filter empty lines, get unique, sort
            ids = sorted(list(set(line.strip() for line in f if line.strip())))
        logging.info(f"Read {len(ids)} unique UniProt IDs from {filename}.")
        return ids
    except Exception as e:
        logging.error(f"Error reading input file {filename}: {e}")
        raise

def load_processed_ids(filename: Path) -> set[str]:
    """Loads the set of already processed UniProt IDs from the TSV output file."""
    processed: set[str] = set()
    if filename.exists() and filename.stat().st_size > 0: # Check size > 0
        logging.info(f"Found existing results file: {filename}. Loading processed IDs...")
        try:
            with open(filename, 'r', newline='', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\\t')
                try:
                    header = next(reader) # Read header
                    if header != ["UniProtKB_AC", "PDB_IDs"]:
                         logging.warning(f"Output file header mismatch. Expected ['UniProtKB_AC', 'PDB_IDs'], got {header}. Attempting to read anyway.")
                except StopIteration:
                    logging.warning(f"Output file {filename} exists but is empty or has no header. Starting fresh.")
                    return set() # File exists but is empty

                # Read the first column (UniProtKB_AC)
                for row in reader:
                    if row and row[0]: # Check if row and first element exist
                        processed.add(row[0])
            logging.info(f"Loaded {len(processed)} previously processed IDs.")
        except Exception as e:
            logging.warning(f"Error reading existing results file {filename}: {e}. Assuming no IDs processed.")
            return set() # Return empty set on error
    else:
        logging.info(f"No existing results file found or file is empty ({filename}).")
    return processed

def query_uniprot_single(accession: str, max_retries: int, retry_delay: int) -> tuple[str, str]:
    """
    Queries the UniProt API for a single accession using GET.

    Args:
        accession: The UniProtKB accession to query.
        max_retries: Max retries for transient errors (e.g., timeout, 5xx).
        retry_delay: Delay between retries in seconds.

    Returns:
        Tuple (accession, result_string). result_string contains comma-separated
        PDB IDs, an empty string if none found, or an error message starting with "Error:".
    """
    retries = 0
    api_url = API_URL_TEMPLATE.format(accession=accession)
    params = {'fields': FIELDS_TO_RETRIEVE, 'format': RESPONSE_FORMAT}
    logging.debug(f"Querying URL: {api_url} with params: {params}")

    while retries <= max_retries:
        try:
            response = requests.get(api_url, params=params, timeout=30) # 30s timeout

            # Handle specific non-retryable client errors first
            if response.status_code == 400:
                 logging.warning(f"  Received 400 Bad Request for {accession}. Invalid AC? Skipping.")
                 return (accession, "Error: 400 Bad Request")
            if response.status_code == 404:
                 logging.info(f"  Received 404 Not Found for {accession} (Accession invalid or not found).")
                 return (accession, "Error: 404 Accession Not Found")
            if response.status_code == 410:
                 logging.info(f"  Received 410 Gone for {accession} (Accession obsolete).")
                 return (accession, "Error: 410 Obsolete Accession")

            # Raise errors for other bad responses (e.g., 5xx server errors, other 4xx)
            response.raise_for_status()

            # --- Successful response (200 OK) ---
            # Use StringIO to treat the response text like a file for csv.reader
            tsv_io = io.StringIO(response.text)
            reader = csv.reader(tsv_io, delimiter='\\t')

            try:
                header = next(reader) # Read header line
                logging.debug(f"  Header received for {accession}: {header}")
                # We expect ['Entry', 'Cross-reference (PDB)'] or just ['PDB'] if only field requested
                # No strict check needed, just try to read the next line

                try:
                    data_line = next(reader) # Read the data line
                    logging.debug(f"  Data line received for {accession}: {data_line}")
                    if len(data_line) >= 2: # Should have at least Entry and PDB field
                        pdb_string = data_line[1].strip()
                        # Clean the PDB string: remove trailing semicolons if present
                        pdb_string = pdb_string.strip(';')
                        # Replace internal semicolons with commas for consistency? Optional.
                        # pdb_string = pdb_string.replace(';', ',')
                        return (accession, pdb_string)
                    elif len(data_line) == 1: # Entry found, but no PDB field returned
                         logging.debug(f"  Data line only contained entry for {accession}, assuming no PDB refs.")
                         return (accession, "") # Return empty string for PDB IDs
                    else: # Unexpected data format
                        logging.warning(f"Unexpected data format (empty line?) for {accession}")
                        return (accession, "Error: Unexpected TSV Format")

                except StopIteration:
                    # No data line found after header -> AC is valid but no PDB entry found
                    logging.debug(f"  No data line found for {accession} after header. No PDB refs.")
                    return (accession, "") # Return empty string

            except StopIteration:
                 # Response was completely empty (no header even)? Should not happen on 200 OK.
                 logging.warning(f"Empty response received despite 200 OK for {accession}")
                 return (accession, "Error: Empty Response (200 OK)")

        except requests.exceptions.Timeout:
            retries += 1
            logging.warning(f"  Request timed out for {accession} (Attempt {retries}/{max_retries+1}). Retrying in {retry_delay}s...")
            if retries <= max_retries: time.sleep(retry_delay)
        except requests.exceptions.RequestException as e: # Includes HTTPError for 5xx etc.
            retries += 1
            logging.warning(f"  Request failed for {accession} (Attempt {retries}/{max_retries+1}): {e}. Retrying in {retry_delay}s...")
            if retries <= max_retries: time.sleep(retry_delay)
        except Exception as e:
             logging.error(f"  Unexpected error during UniProt query for {accession}: {e}")
             logging.exception("Traceback:") # Log traceback for unexpected errors
             return (accession, f"Error: Unexpected {type(e).__name__}")

    # If loop finishes without success
    logging.error(f"Query failed for {accession} after {max_retries+1} attempts.")
    return (accession, "Error: Max Retries Exceeded")


def save_single_result(result_tuple: tuple[str, str] | None, filename: Path, write_header: bool):
    """Appends a single result tuple (ac, pdb_string) to the TSV output file."""
    if not result_tuple:
        logging.warning("Attempted to save an empty result tuple. Skipping.")
        return False

    accession, pdb_string = result_tuple
    try:
        # Open in append mode, ensure newline handling is correct
        with open(filename, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\\t', lineterminator='\\n', quoting=csv.QUOTE_MINIMAL)
            if write_header:
                writer.writerow(["UniProtKB_AC", "PDB_IDs"])
            writer.writerow([accession, pdb_string])
        return True
    except IOError as e:
        logging.error(f"Could not write result for {accession} to file {filename}: {e}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error saving result for {accession}: {e}")
        return False


def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Query UniProtKB ACs for PDB cross-references.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_file", required=True, type=Path,
                        help="Input file containing UniProtKB ACs (one per line).")
    parser.add_argument("-o", "--output_file", required=True, type=Path,
                        help="Output file to save results (TSV format, resumable).")
    parser.add_argument("--delay", type=float, default=0.5,
                        help="Delay between API requests in seconds.")
    parser.add_argument("--max_retries", type=int, default=2, # Increased default slightly
                        help="Max retries for a failed request (timeouts, server errors).")
    parser.add_argument("--retry_delay", type=int, default=5,
                        help="Delay between retries in seconds.")
    parser.add_argument("--progress", type=int, default=100,
                        help="Print progress update every N IDs.")

    return parser.parse_args()

# --- Main Execution ---
if __name__ == "__main__":
    args = parse_arguments()

    # Check dependencies
    if not shutil.which("python3"): # Basic check if python is running
        logging.critical("Python interpreter not found?")
        sys.exit(1)
    try:
        import requests
    except ImportError:
        logging.critical("The 'requests' library is required. Please install it (`pip install requests`).")
        sys.exit(1)


    try:
        all_target_ids = read_uniprot_ids(args.input_file)
        if not all_target_ids:
            logging.info("Input file is empty or no valid IDs found. Exiting.")
            sys.exit(0)

        processed_ids_set = load_processed_ids(args.output_file)
        ids_to_query_now = [pid for pid in all_target_ids if pid not in processed_ids_set]
        total_remaining = len(ids_to_query_now)

        if not ids_to_query_now:
            logging.info("All UniProt IDs from the input file have already been processed according to the results file.")
            sys.exit(0)

        logging.info(f"Starting query for {total_remaining:,} remaining UniProt IDs (one by one using GET)...")
        # Determine if header needs to be written
        is_new_file = not args.output_file.exists() or args.output_file.stat().st_size == 0

        session_processed_count = 0
        session_saved_count = 0
        session_error_count = 0

        # --- Main Loop ---
        for i, current_id in enumerate(ids_to_query_now):
            # Print progress before making the request
            if (i) % args.progress == 0 and i > 0:
                 logging.info(f"\n--- Progress: Attempted {session_processed_count} queries this session ({i}/{total_remaining} remaining IDs) ---")

            logging.info(f"Processing ID {i+1}/{total_remaining}: {current_id}...")
            result = query_uniprot_single(current_id, args.max_retries, args.retry_delay)

            # Process result
            if result:
                acc, res_string = result
                if "Error:" in res_string:
                    logging.warning(f"  Result for {acc}: {res_string}") # Log errors as warnings
                    session_error_count += 1
                elif not res_string:
                    logging.info(f"  Result for {acc}: No PDB Found")
                else:
                    logging.info(f"  Result for {acc}: Found PDB: {res_string[:100]}{'...' if len(res_string)>100 else ''}")

                # Save result regardless of error status (to mark as processed)
                saved_ok = save_single_result(result, args.output_file, is_new_file)
                if saved_ok:
                    session_saved_count += 1
                    is_new_file = False # Header is written only once
                else:
                    logging.error(f"  Failed to save result for {current_id}")
                    session_error_count += 1 # Count save failure as an error
            else:
                 # Should not happen if query_uniprot_single always returns a tuple
                 logging.error(f"  No result tuple obtained for {current_id}. This indicates an internal script error.")
                 session_error_count += 1

            session_processed_count += 1

            # Delay between requests (only if more IDs remain)
            if i + 1 < total_remaining:
                time.sleep(args.delay)
        # --- End Loop ---

        # Final summary
        final_processed_in_file = len(load_processed_ids(args.output_file)) # Re-check file
        logging.info(f"\n--- Processing Complete ---")
        logging.info(f"Attempted queries for {session_processed_count} IDs in this session.")
        logging.info(f"Total unique IDs processed and saved in output file: {final_processed_in_file}/{len(all_target_ids)}")
        logging.info(f"Total result lines saved/appended in this session: {session_saved_count}")
        logging.info(f"Total errors encountered (query or save): {session_error_count}")
        logging.info(f"Output file: {args.output_file}")

        if session_error_count > 0:
             logging.warning("Script finished with errors. Check log messages.")
             sys.exit(1) # Exit with error code
        else:
             logging.info("Script finished successfully.")
             sys.exit(0)

    except FileNotFoundError as e:
        logging.error(e)
        sys.exit(1)
    except Exception as e:
        logging.critical(f"An unexpected critical error occurred during execution: {e}")
        logging.exception("Traceback:")
        sys.exit(1)

