import pandas as pd
import requests
import json
import time
import sys
import argparse
import logging
from pathlib import Path
from typing import Union # Import Union for older Python compatibility
from tqdm.auto import tqdm # Optional: for progress bar (pip install tqdm)

# --- Setup Logging ---
# Configured in main() function via argparse

# --- Constants ---
AFDB_JSON_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{identifier}-F1-confidence_v4.json"
DEFAULT_REQUEST_DELAY = 0.1

def fetch_and_calculate_plddt(protein_id_debug: str, identifier: str, afdb_url_template: str) -> Union[float, None]:
    """
    Fetches the AFDB JSON (from confidence endpoint) for a given identifier 
    and calculates the average pLDDT (using 'confidenceScore' key).

    Args:
        protein_id_debug: The internal ProteinID for debugging messages.
        identifier: The UniProt AC or UPI.
        afdb_url_template: URL template for fetching AFDB data.

    Returns:
        The average pLDDT score as a float, or None if fetching/calculation fails.
    """
    calculated_avg_plddt = None 

    if pd.isna(identifier):
        logging.debug(f"[{protein_id_debug}]: Identifier is NaN. Skipping.")
        logging.info(f"[{protein_id_debug} | NaN_Identifier]: Failed to retrieve Avg pLDDT. Reason: Identifier was NaN.")
        return None
    
    identifier_str = str(identifier).strip()

    if not identifier_str:
        logging.debug(f"[{protein_id_debug}]: Identifier is empty after stripping. Skipping.")
        logging.info(f"[{protein_id_debug} | Empty_Identifier]: Failed to retrieve Avg pLDDT. Reason: Identifier was empty.")
        return None
    
    url = afdb_url_template.format(identifier=identifier_str)
    logging.debug(f"[{protein_id_debug} | {identifier_str}]: Processing AFDB ID. Fetching URL: {url}")

    try:
        response = requests.get(url, timeout=30)
        logging.debug(f"[{protein_id_debug} | {identifier_str}]: Response Status Code: {response.status_code}")
        response.raise_for_status() 

        data = response.json()
        
        logging.debug(f"[{protein_id_debug} | {identifier_str}]: Type of 'data' object: {type(data)}")
        logging.debug(f"[{protein_id_debug} | {identifier_str}]: Raw 'data' (first 500 chars): {str(data)[:500]}")
        if isinstance(data, dict):
            logging.debug(f"[{protein_id_debug} | {identifier_str}]: Keys in 'data' dict: {list(data.keys())}")
        elif isinstance(data, list) and data: 
            logging.debug(f"[{protein_id_debug} | {identifier_str}]: 'data' is a list with {len(data)} element(s). Type of 'data[0]': {type(data[0])}")
            if isinstance(data[0], dict):
                logging.debug(f"[{protein_id_debug} | {identifier_str}]: Keys in 'data[0]' dict: {list(data[0].keys())}")

        confidence_scores = None
        if isinstance(data, dict): 
            confidence_scores = data.get('confidenceScore') 
        elif isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict): 
            confidence_scores = data[0].get('confidenceScore')
        else:
            logging.warning(f"[{protein_id_debug} | {identifier_str}]: JSON 'data' object is not a dict or a list of dicts as expected.")

        logging.debug(f"[{protein_id_debug} | {identifier_str}]: Extracted 'confidence_scores' (type: {type(confidence_scores)}): {str(confidence_scores)[:200] if confidence_scores is not None else 'None'}")

        if confidence_scores and isinstance(confidence_scores, list) and len(confidence_scores) > 0:
            if all(isinstance(score, (int, float)) for score in confidence_scores):
                avg_confidence = sum(confidence_scores) / len(confidence_scores)
                calculated_avg_plddt = round(avg_confidence, 2)
                logging.debug(f"[{protein_id_debug} | {identifier_str}]: Successfully calculated Avg Confidence Score (pLDDT): {calculated_avg_plddt}")
            else:
                logging.warning(f"[{protein_id_debug} | {identifier_str}]: 'confidenceScore' list contains non-numeric values. Scores: {str(confidence_scores[:10])}")
        else:
            logging.debug(f"[{protein_id_debug} | {identifier_str}]: Final check: 'confidenceScore' is None, not a list, or an empty list.")

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            logging.info(f"[{protein_id_debug} | {identifier_str}]: AFDB confidence entry not found (404).")
        else:
            logging.error(f"[{protein_id_debug} | {identifier_str}]: HTTP Error {e.response.status_code} - {e.response.text[:200]}.")
    except requests.exceptions.RequestException as e:
        logging.error(f"[{protein_id_debug} | {identifier_str}]: Request failed: {e}.")
    except json.JSONDecodeError as e_json:
        logging.error(f"[{protein_id_debug} | {identifier_str}]: Could not decode JSON. Error: {e_json}.")
    except Exception as e:
        logging.error(f"[{protein_id_debug} | {identifier_str}]: An unexpected error occurred: {e}.")
        logging.exception("Traceback for unexpected error:")
    
    if calculated_avg_plddt is not None:
        logging.info(f"[{protein_id_debug} | {identifier_str}]: Avg Confidence Score (pLDDT) {calculated_avg_plddt} successfully retrieved from {url}")
    else:
        logging.info(f"[{protein_id_debug} | {identifier_str}]: Failed to retrieve Avg Confidence Score (pLDDT) from {url}.")
        
    return calculated_avg_plddt

def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Fetch pLDDT scores from AlphaFold DB for identifiers in a CSV file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_csv", required=True, type=Path,
                        help="Input CSV file containing protein identifiers.")
    parser.add_argument("-o", "--output_csv", required=True, type=Path,
                        help="Output CSV file to save pLDDT scores.")
    parser.add_argument("--protein_id_col", default="ProteinID",
                        help="Column name in input CSV for your internal protein ID.")
    parser.add_argument("--afdb_id_col", required=True,
                        help="Column name in input CSV for the UniProt AC or UPI to query AFDB.")
    parser.add_argument("--delay", type=float, default=DEFAULT_REQUEST_DELAY,
                        help="Delay in seconds between API requests.")
    parser.add_argument("--afdb_url_template", default=AFDB_JSON_URL_TEMPLATE,
                        help="URL template for AlphaFold DB JSON files.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Optional: Process only the first N rows for testing.")
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Set the logging level.")
    return parser.parse_args()

# --- Main Script ---
def main():
    args = parse_arguments()

    # Setup logging based on command line argument
    numeric_log_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError(f"Invalid log level: {args.log_level}")
    logging.basicConfig(level=numeric_log_level,
                        format='%(asctime)s - %(levelname)s - [%(funcName)s] - %(message)s',
                        stream=sys.stdout)

    logging.info(f"Starting pLDDT fetch process...")
    logging.info(f"Input file: {args.input_csv}")
    logging.info(f"Output file: {args.output_csv}")
    logging.info(f"Protein ID column: {args.protein_id_col}")
    logging.info(f"AFDB ID column: {args.afdb_id_col}")
    logging.info(f"Request delay: {args.delay}s")

    results = []

    try:
        df_input = pd.read_csv(args.input_csv, low_memory=False)
        logging.info(f"Read {len(df_input)} rows from {args.input_csv}.")

        if args.protein_id_col not in df_input.columns:
             logging.error(f"Error: Column '{args.protein_id_col}' not found in {args.input_csv}. Please check configuration.")
             sys.exit(1)
        if args.afdb_id_col not in df_input.columns:
             logging.error(f"Error: Column '{args.afdb_id_col}' not found in {args.input_csv}. Please check configuration.")
             sys.exit(1)

        logging.info(f"\nFirst 5 AFDB Identifiers to be processed (from column '{args.afdb_id_col}'):")
        for i, af_id in enumerate(df_input[args.afdb_id_col].head()):
            logging.info(f"  {i+1}: '{af_id}' (type: {type(af_id)})")
        
        df_to_process = df_input
        if args.limit is not None and args.limit > 0:
            logging.info(f"Processing only the first {args.limit} rows due to --limit option.")
            df_to_process = df_input.head(args.limit)

        logging.info(f"\nFetching pLDDT scores from AlphaFold DB (Processing {len(df_to_process)} entries)...")

        # Ensure output directory exists
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)

        # Use tqdm for progress bar if available
        iterable_rows = df_to_process.iterrows()
        if 'tqdm' in sys.modules:
             iterable_rows = tqdm(df_to_process.iterrows(), total=len(df_to_process), desc="Processing Proteins")

        for index, row in iterable_rows:
            protein_id = row[args.protein_id_col]
            afdb_id_original = row[args.afdb_id_col]

            avg_plddt = fetch_and_calculate_plddt(str(protein_id), afdb_id_original, args.afdb_url_template)

            results.append({
                args.protein_id_col: protein_id, # Use the specified column name
                'Avg_pLDDT': avg_plddt
            })
            time.sleep(args.delay)

    except FileNotFoundError:
        logging.error(f"Error: Input file not found at '{args.input_csv}'")
        sys.exit(1)
    except Exception as e:
        logging.critical(f"An unexpected error occurred during processing: {e}")
        logging.exception("Traceback for critical error:")
        sys.exit(1)

    if results:
        df_results = pd.DataFrame(results)
        successful_fetches = df_results['Avg_pLDDT'].notna().sum()
        logging.info(f"\nSuccessfully fetched {successful_fetches} out of {len(df_to_process)} pLDDT scores.")
        
        try:
            df_results.to_csv(args.output_csv, index=False, na_rep='NA')
            logging.info(f"Successfully saved results to {args.output_csv}")
        except Exception as e:
            logging.error(f"Error writing output CSV file: {e}")
    else:
        logging.info("No results were generated to save.")

    logging.info("Script finished.")

if __name__ == "__main__":
    main()
