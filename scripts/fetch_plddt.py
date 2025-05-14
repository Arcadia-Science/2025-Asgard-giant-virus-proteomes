import pandas as pd
import requests
import json
import time
import sys
from typing import Union # Import Union for older Python compatibility
from tqdm.auto import tqdm # Optional: for progress bar (pip install tqdm)

# --- Configuration ---
# Input CSV file containing ProteinIDs and their corresponding AFDB identifiers (UniProt ACs/UPIs)
INPUT_CSV = "PDB_AFDB_screening/afdb_processing_output/afdb_found_uniprot_acs_or_upi.csv" 

# Output CSV file path
OUTPUT_CSV = "protein_avg_plddt_scores_confidence_endpoint_corrected_key.csv" # New output name

# Column name in INPUT_CSV that holds your internal Protein ID
PROTEIN_ID_COLUMN = "Uniprot_AC" 

# Column name in INPUT_CSV that holds the UniProt AC or UPI identifier found in AFDB
AFDB_IDENTIFIER_COLUMN = "Uniprot_AC" # Adjust if necessary

# Base URL for AlphaFold DB JSON files - NOW USING THE CONFIDENCE ENDPOINT
AFDB_JSON_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{identifier}-F1-confidence_v4.json"

# Delay between requests in seconds (to be polite to the server)
REQUEST_DELAY = 0.1 
# --- End Configuration ---

def fetch_and_calculate_plddt(protein_id_debug: str, identifier: str) -> Union[float, None]:
    """
    Fetches the AFDB JSON (from confidence endpoint) for a given identifier 
    and calculates the average pLDDT (now using 'confidenceScore' key). 
    Includes enhanced debugging print statements and a summary INFO line for each attempt.

    Args:
        protein_id_debug: The internal ProteinID for debugging messages.
        identifier: The UniProt AC or UPI.

    Returns:
        The average pLDDT score as a float, or None if fetching/calculation fails.
    """
    calculated_avg_plddt = None 

    if pd.isna(identifier):
        print(f"DEBUG [{protein_id_debug}]: Identifier is NaN. Skipping.")
        print(f"INFO [{protein_id_debug} | NaN_Identifier]: Failed to retrieve Avg pLDDT. Reason: Identifier was NaN.")
        return None
    
    identifier_str = str(identifier).strip()

    if not identifier_str:
        print(f"DEBUG [{protein_id_debug}]: Identifier is empty after stripping. Skipping.")
        print(f"INFO [{protein_id_debug} | Empty_Identifier]: Failed to retrieve Avg pLDDT. Reason: Identifier was empty.")
        return None
    
    url = AFDB_JSON_URL_TEMPLATE.format(identifier=identifier_str)
    print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Processing AFDB ID.")
    print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Fetching URL: {url}")

    try:
        response = requests.get(url, timeout=30)
        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Response Status Code: {response.status_code}")
        response.raise_for_status() 

        data = response.json()
        
        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Type of 'data' object: {type(data)}")
        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Raw 'data' (first 500 chars): {str(data)[:500]}")
        if isinstance(data, dict):
            print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Keys in 'data' dict: {list(data.keys())}")
        elif isinstance(data, list) and data: 
            print(f"DEBUG [{protein_id_debug} | {identifier_str}]: 'data' is a list with {len(data)} element(s). Type of 'data[0]': {type(data[0])}")
            if isinstance(data[0], dict):
                print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Keys in 'data[0]' dict: {list(data[0].keys())}")

        confidence_scores = None # Changed variable name for clarity
        if isinstance(data, dict): 
            # CORRECTED KEY: Using 'confidenceScore' instead of 'plddt'
            confidence_scores = data.get('confidenceScore') 
        elif isinstance(data, list) and len(data) > 0 and isinstance(data[0], dict): 
            # CORRECTED KEY: Using 'confidenceScore' instead of 'plddt'
            confidence_scores = data[0].get('confidenceScore')
        else:
            print(f"DEBUG [{protein_id_debug} | {identifier_str}]: JSON 'data' object is not a dict or a list of dicts as expected.")

        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Extracted 'confidence_scores' (type: {type(confidence_scores)}): {str(confidence_scores)[:200] if confidence_scores is not None else 'None'}")

        if confidence_scores and isinstance(confidence_scores, list) and len(confidence_scores) > 0:
            if all(isinstance(score, (int, float)) for score in confidence_scores):
                avg_confidence = sum(confidence_scores) / len(confidence_scores)
                calculated_avg_plddt = round(avg_confidence, 2)
                print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Successfully calculated Avg Confidence Score (pLDDT): {calculated_avg_plddt}")
            else:
                print(f"DEBUG [{protein_id_debug} | {identifier_str}]: 'confidenceScore' list contains non-numeric values. Scores: {confidence_scores[:10]}")
        else:
            print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Final check: 'confidenceScore' is None, not a list, or an empty list.")

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            print(f"DEBUG [{protein_id_debug} | {identifier_str}]: AFDB confidence entry not found (404).")
        else:
            print(f"DEBUG [{protein_id_debug} | {identifier_str}]: HTTP Error {e.response.status_code} - {e.response.text[:200]}.")
    except requests.exceptions.RequestException as e:
        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Request failed: {e}.")
    except json.JSONDecodeError as e_json:
        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: Could not decode JSON. Error: {e_json}.")
    except Exception as e:
        print(f"DEBUG [{protein_id_debug} | {identifier_str}]: An unexpected error occurred: {e}.")
        import traceback
        traceback.print_exc()
    
    if calculated_avg_plddt is not None:
        print(f"INFO [{protein_id_debug} | {identifier_str}]: Avg Confidence Score (pLDDT) {calculated_avg_plddt} successfully retrieved from {url}")
    else:
        print(f"INFO [{protein_id_debug} | {identifier_str}]: Failed to retrieve Avg Confidence Score (pLDDT) from {url}. (See DEBUG logs for details)")
        
    return calculated_avg_plddt

# --- Main Script ---
if __name__ == "__main__":
    print(f"Starting pLDDT fetch process (CONFIDENCE endpoint, CORRECTED JSON KEY)...") 
    print(f"Input file: {INPUT_CSV}")
    print(f"Output file: {OUTPUT_CSV}")

    results = []

    try:
        df_input = pd.read_csv(INPUT_CSV, low_memory=False)
        print(f"Read {len(df_input)} rows from {INPUT_CSV}.")

        if PROTEIN_ID_COLUMN not in df_input.columns:
             print(f"Error: Column '{PROTEIN_ID_COLUMN}' not found in {INPUT_CSV}. Please check configuration.")
             sys.exit(1)
        if AFDB_IDENTIFIER_COLUMN not in df_input.columns:
             print(f"Error: Column '{AFDB_IDENTIFIER_COLUMN}' not found in {INPUT_CSV}. Please check configuration.")
             sys.exit(1)

        print(f"\nFirst 5 AFDB Identifiers to be processed (from column '{AFDB_IDENTIFIER_COLUMN}'):")
        for i, af_id in enumerate(df_input[AFDB_IDENTIFIER_COLUMN].head()):
            print(f"  {i+1}: '{af_id}' (type: {type(af_id)})")
        print("-" * 30)
        
        # Optional: Process only a small subset for initial debugging
        # df_to_process = df_input.head(5) 
        df_to_process = df_input

        print(f"\nFetching pLDDT scores from AlphaFold DB (Processing {len(df_to_process)} entries)...")
        for index, row in tqdm(df_to_process.iterrows(), total=len(df_to_process), desc="Processing Proteins"):
            protein_id = row[PROTEIN_ID_COLUMN]
            afdb_id_original = row[AFDB_IDENTIFIER_COLUMN]

            avg_plddt = fetch_and_calculate_plddt(str(protein_id), afdb_id_original) 

            results.append({
                PROTEIN_ID_COLUMN: protein_id,
                'Avg_pLDDT': avg_plddt # Still calling the output column Avg_pLDDT for consistency
            })
            time.sleep(REQUEST_DELAY)

    except FileNotFoundError:
        print(f"Error: Input file not found at '{INPUT_CSV}'")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during processing: {e}")
        import traceback
        traceback.print_exc() 
        sys.exit(1)

    if results:
        df_results = pd.DataFrame(results)
        successful_fetches = df_results['Avg_pLDDT'].notna().sum()
        print(f"\nSuccessfully fetched {successful_fetches} out of {len(df_to_process)} pLDDT scores.")
        
        try:
            df_results.to_csv(OUTPUT_CSV, index=False, na_rep='NA')
            print(f"Successfully saved results to {OUTPUT_CSV}")
        except Exception as e:
            print(f"Error writing output CSV file: {e}")
    else:
        print("No results were generated.")

    print("Script finished.")
