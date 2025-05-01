#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Searches the UniParc REST API (/search endpoint) for a list of protein IDs
(e.g., RefSeq, GenBank CDS accessions) to find corresponding UniParc
Identifiers (UPIs).

This is typically used for IDs that failed to map directly to UniProtKB.
"""

import requests
import time
import sys
from pathlib import Path
import csv
import traceback
import argparse
import datetime
import logging
import json

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# --- Constants ---
UNIPARC_SEARCH_URL = "https://rest.uniprot.org/uniparc/search"
DEFAULT_SLEEP_INTERVAL = 0.5 # Seconds between requests

# --- Helper Function ---

def search_uniparc_for_id(protein_id: str) -> tuple[str | None, str | None]:
    """
    Searches UniParc using the /search endpoint for a given protein ID.

    Args:
        protein_id: The protein ID (e.g., RefSeq, GenBank accession).

    Returns:
        Tuple (upi, error_message).
        - upi (str | None): The UniParc Identifier (UPI) string if found, None otherwise.
        - error_message (str | None): None if successful (including 'Not Found'),
          or a string describing the error otherwise.
    """
    query_params = {
        "query": f"{protein_id}", # Search for the ID itself
        "format": "tsv",
        "fields": "upi", # Only retrieve the UPI field
        "size": 1 # We only need the first hit if found
    }
    url = UNIPARC_SEARCH_URL
    logging.debug(f"Searching UniParc for ID: {protein_id} with query: {query_params}")

    try:
        response = requests.get(url, params=query_params, timeout=30) # 30s timeout
        response.raise_for_status() # Check for HTTP errors (4xx, 5xx)

        # Response is expected to be TSV
        content = response.text.strip()
        logging.debug(f"  Response text for {protein_id}: '{content[:100]}...'")

        # Split lines and parse TSV content
        lines = content.splitlines()

        # Check if we got results (more than just a header line)
        if len(lines) > 1:
            # Assume the first column on the second line is the UPI
            try:
                # Use csv reader for safer TSV parsing
                reader = csv.reader([lines[1]], delimiter='\\t') # Read only the second line
                data_row = next(reader)
                if data_row: # Check if row is not empty
                    possible_upi = data_row[0].strip()
                    if possible_upi.startswith("UPI"):
                        logging.debug(f"  Found UPI: {possible_upi} for {protein_id}")
                        return possible_upi, None # Success
                    else:
                        error_msg = f"Unexpected content in first column of data line: '{possible_upi}' for ID {protein_id}. Full line: '{lines[1]}'"
                        logging.warning(error_msg)
                        return None, error_msg # Treat as error for logging
                else:
                     error_msg = f"Empty data row parsed for ID {protein_id}. Line content: '{lines[1]}'"
                     logging.warning(error_msg)
                     return None, error_msg # Treat as error

            except (IndexError, StopIteration) as parse_e:
                 # Handle cases where the second line might be empty or TSV split failed
                 error_msg = f"Could not parse data line for ID {protein_id}. Line content: '{lines[1] if len(lines) > 1 else 'N/A'}'. Error: {parse_e}"
                 logging.warning(error_msg)
                 return None, error_msg # Treat as error

        elif len(lines) == 1:
             # Only a header line returned - means not found
             logging.debug(f"  Only header line returned for {protein_id}: '{lines[0]}'. ID not found.")
             return None, None # Treat as not found
        else: # Empty content
             logging.debug(f"  Empty response content for {protein_id}. ID not found.")
             return None, None # Treat as not found

    except requests.exceptions.Timeout as e:
        error_msg = f"Timeout searching UniParc for {protein_id}: {e}"
        logging.error(error_msg)
        return None, error_msg
    except requests.exceptions.HTTPError as e:
        status_code = e.response.status_code
        response_text = e.response.text[:200]
        error_msg = f"HTTP Error {status_code} searching UniParc for {protein_id}. Response: {response_text}"
        logging.error(error_msg)
        return None, error_msg # Treat all HTTP errors as lookup failures
    except requests.exceptions.RequestException as e:
        error_msg = f"Network error searching UniParc for {protein_id}: {e}"
        logging.error(error_msg)
        return None, error_msg
    except Exception as e:
        error_msg = f"Unexpected error searching UniParc for {protein_id}: {e}"
        logging.exception(error_msg) # Log traceback for unexpected errors
        return None, error_msg

# --- Argument Parser ---
def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Search UniParc for protein IDs using the /search endpoint.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_file", required=True, type=Path,
                        help="Input file containing protein IDs (one per line).")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Directory for output files (mapping TSV, not found list, error log).")
    parser.add_argument("--output_prefix", default=None,
                        help="Optional prefix for output files. If None, uses input file stem.")
    parser.add_argument("--sleep", type=float, default=DEFAULT_SLEEP_INTERVAL,
                        help="Seconds to wait between API requests.")

    return parser.parse_args()

# --- Main Script Logic ---
def main():
    """Main function to search UniParc by ID using the search endpoint."""
    args = parse_arguments()

    logging.info("--- Starting UniParc Search by ID ---")
    logging.info(f"Reading target IDs from: {args.input_file}")

    # --- Step 1: Check and Load Target IDs ---
    target_ids = []
    try:
        with open(args.input_file, 'r') as f:
            target_ids = [line.strip() for line in f if line.strip()]
        logging.info(f"Loaded {len(target_ids):,} IDs to search in UniParc.")
        if not target_ids:
            logging.error("Error: No IDs found in the input file.")
            sys.exit(1)
    except FileNotFoundError:
        logging.error(f"Error: Input IDs file not found: {args.input_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading IDs file {args.input_file}: {e}")
        sys.exit(1)

    # --- Step 2: Prepare Output Files ---
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create specific output filenames
    if args.output_prefix:
        output_base = args.output_prefix
    else:
        output_base = args.input_file.stem # Use the stem of the input file
    output_mapping_file = output_dir / f"{output_base}_uniparc_mapping.tsv"
    output_not_found_file = output_dir / f"{output_base}_uniparc_notfound.txt"
    output_error_file = output_dir / f"{output_base}_uniparc_errors.json" # Log errors to JSON

    logging.info(f"Output mapping file: {output_mapping_file}")
    logging.info(f"Not found IDs file: {output_not_found_file}")
    logging.info(f"Error log file: {output_error_file}")

    # --- Step 3: Iterate and Search UniParc ---
    successful_mappings = []
    not_found_ids = []
    error_details = [] # Store dicts of {'id': id, 'error': msg}
    total_ids = len(target_ids)

    logging.info(f"\n--- Searching UniParc for {total_ids:,} IDs ---")
    try:
        # Open output file in write mode, with newline='' for csv module
        with open(output_mapping_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile, delimiter='\\t', lineterminator='\\n')
            writer.writerow(["Input_ID", "UniParc_ID"]) # Write header

            # Loop through each ID from the input file
            for i, input_id in enumerate(target_ids):
                # Log progress periodically
                if (i + 1) % 100 == 0 or (i + 1) == total_ids:
                     logging.info(f"Processing {i+1}/{total_ids}: {input_id}...")

                # Call the lookup function
                upi, error_msg = search_uniparc_for_id(input_id)

                # Process the result
                if upi:
                    # Success: Found a UPI
                    logging.debug(f"  Found UPI: {upi} for {input_id}")
                    writer.writerow([input_id, upi])
                    outfile.flush() # Write mapping immediately
                    successful_mappings.append(input_id)
                elif error_msg:
                    # Error occurred during lookup
                    logging.warning(f"  Error searching for {input_id}: {error_msg}")
                    error_details.append({"id": input_id, "error": error_msg})
                    # Decide if errors should also be added to not_found list
                    # not_found_ids.append(input_id)
                else:
                    # No UPI found (API returned no results, 404, etc.)
                    logging.debug(f"  ID not found in UniParc search: {input_id}")
                    not_found_ids.append(input_id)

                # Rate limiting between API calls
                time.sleep(args.sleep)

    except IOError as e:
         logging.critical(f"Error writing to output mapping file {output_mapping_file}: {e}")
         sys.exit(1) # Exit if we can't write output
    except Exception as e:
        logging.critical(f"An critical error occurred during the search process: {e}")
        logging.exception("Traceback:")
        logging.critical("The search process may be incomplete.")

    # --- Step 4: Save Not Found and Error IDs ---
    if not_found_ids:
        logging.info(f"\nSaving {len(not_found_ids)} IDs not found in UniParc search to: {output_not_found_file}")
        try:
            with open(output_not_found_file, 'w') as f_not_found:
                for prot_id in sorted(not_found_ids):
                    f_not_found.write(f"{prot_id}\\n")
        except Exception as e:
            logging.error(f"Error writing not found IDs file: {e}")

    if error_details:
        logging.warning(f"\nSaving details for {len(error_details)} IDs that encountered errors to: {output_error_file}")
        try:
            with open(output_error_file, 'w') as f_err:
                json.dump(error_details, f_err, indent=2)
        except Exception as e:
             logging.error(f"Error writing error log file: {e}")

    # --- Final Summary ---
    logging.info("\n--- UniParc Search by ID Finished ---")
    logging.info(f"Total IDs processed: {total_ids:,}")
    logging.info(f"Successful searches finding a UPI: {len(successful_mappings):,}")
    logging.info(f"IDs not found via UniParc search: {len(not_found_ids):,}")
    logging.info(f"IDs resulting in search errors: {len(error_details):,}")
    logging.info(f"Mapping results saved to: {output_mapping_file}")
    if not_found_ids:
        logging.info(f"IDs not found are listed in: {output_not_found_file}")
    if error_details:
        logging.warning(f"Error details saved to: {output_error_file}")

    if error_details:
        sys.exit(1) # Exit with error code if errors occurred
    else:
        sys.exit(0)


# --- Script Entry Point ---
if __name__ == "__main__":
    main()
