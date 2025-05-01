#!/usr/bin/env python3

"""
Selects potential outgroup sequences (e.g., TACK, Euryarchaeota) for multiple
Orthologous Groups (OGs) based on DIAMOND search results.

Reads DIAMOND output files (one per OG), classifies hits based on keywords found
within the sequence ID (sseqid), filters hits based on E-value and coverage,
and selects a defined number of top-scoring hits, prioritizing diversity
(e.g., one TACK, one Eury). Saves the selected outgroup IDs to a summary CSV file.
"""

import pandas as pd
import os
import re
import argparse
import glob
from collections import defaultdict
import sys
import logging
from pathlib import Path

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - [%(funcName)s] %(message)s',
                    stream=sys.stdout)

# --- Default Configuration (Can be overridden or loaded from files) ---
# Filtering Thresholds
DEFAULT_EVALUE_THRESHOLD = 1e-10
DEFAULT_MIN_COVERAGE = 0.50 # Query and Subject coverage
DEFAULT_MAX_OUTGROUPS_PER_OG = 3

# --- Default Keyword Lists (ADAPT THESE FOR YOUR OUTGROUPS) ---
# Using broader taxonomic levels might be more robust than just genera.
DEFAULT_TACK_KEYWORDS = [
    'Thermoproteota', 'Crenarchaeota', 'Nitrososphaerota', 'Thaumarchaeota',
    'Korarchaeota', 'Aigarchaeota', 'Bathyarchaeota',
    'Thermoprotei', 'Nitrososphaeria', 'Korarchaeia', 'Sulfolobus', 'Thermoproteus',
    'Pyrobaculum', 'Desulfurococcus', 'Ignicoccus', 'Nitrosopumilus', 'Nitrososphaera',
    'Nitrosotalea', 'Nitrosocaldus', 'Korarchaeum', 'Caldiarchaeum', 'Cenarchaeum',
    'Thermosphaera'
]
DEFAULT_EURY_KEYWORDS = [
    'Euryarchaeota', 'Methanobacteria', 'Methanococci', 'Methanomicrobia', 'Halobacteria',
    'Thermococci', 'Archaeoglobi', 'Thermoplasmata', 'Methanopyri', 'Methanobacterium',
    'Methanobrevibacter', 'Methanothermobacter', 'Methanococcus', 'Methanothermococcus',
    'Methanomicrobium', 'Methanoculleus', 'Methanospirillum', 'Methanosarcina',
    'Methanosaeta', 'Methanothrix', 'Halobacterium', 'Halococcus', 'Haloarcula',
    'Halorubrum', 'Haloquadratum', 'Thermococcus', 'Pyrococcus', 'Archaeoglobus',
    'Ferroglobus', 'Thermoplasma', 'Picrophilus', 'Ferroplasma', 'Aciduliprofundum',
    'Methanopyrus', 'Methanosalsum', 'Methanoregulaceae', 'Methanocorpusculum',
    'Methanodesulfokora', 'Thermoplasmatales', 'Methanocellales', 'Methanobacteriota',
    'Methanosphaera', 'Methanolinea', 'Methanomassiliicoccales', 'Natrialbaceae',
    'Natronococcus', 'Halomicrobium', 'Natronoarchaeum', 'Halomarina', 'Halolamina',
    'Salarchaeum', 'Halosegnis', 'Natrarchaeobius', 'Halorussus', 'Halarchaeum',
    'Methanococcoides', 'Methanoplanus', 'Syntrophoarchaeum', 'Methanoperedens',
    'Methanomarinus', 'ANME-1', 'ANME-2', 'Methanotrichaceae', 'Methanophagales',
    'Alkanophagales', 'Nitrosopumilaceae', 'Hadesarchaea', 'Hydrothermarchaeota',
    'Marine_Group_II', 'Marine_Group_III'
]

# --- Functions ---

def get_outgroup_classification(sseqid: str, tack_keywords: set[str], eury_keywords: set[str]) -> str:
    """
    Classifies a sequence ID as TACK, Euryarchaeota, or Unknown based on keywords.
    Searches the ENTIRE sseqid string for keywords (case-insensitive).

    Args:
        sseqid: The sequence ID string from DIAMOND output.
        tack_keywords: Set of keywords indicating TACK group.
        eury_keywords: Set of keywords indicating Euryarchaeota group.

    Returns:
        'TACK', 'Euryarchaeota', or 'Unknown'.
    """
    try:
        sseqid_str = str(sseqid)
        # Use pre-compiled regex patterns for slight efficiency gain if running on huge datasets
        # For simplicity here, use direct search with word boundaries/delimiters

        # Check TACK keywords first
        for keyword in tack_keywords:
            # Regex: Keyword surrounded by start/end, word boundary, or common delimiters (_, |)
            pattern = r'(?:^|\\b|_|\\|)' + re.escape(keyword) + r'(?:$|\\b|_|\\|)'
            if re.search(pattern, sseqid_str, re.IGNORECASE):
                return 'TACK'

        # Check Euryarchaeota keywords if not TACK
        for keyword in eury_keywords:
            pattern = r'(?:^|\\b|_|\\|)' + re.escape(keyword) + r'(?:$|\\b|_|\\|)'
            if re.search(pattern, sseqid_str, re.IGNORECASE):
                return 'Euryarchaeota'

        # Return Unknown if no keywords matched
        logging.debug(f"Could not classify sseqid: {sseqid_str}")
        return 'Unknown'

    except Exception as e:
        logging.warning(f"Error parsing taxonomy for sseqid '{sseqid}': {e}")
        return 'Unknown'


def process_diamond_file(filepath: Path, og_id: str, tack_kw: set[str], eury_kw: set[str],
                         evalue_thr: float, min_cov: float, max_outgroups: int) -> list[str]:
    """ Reads, filters, and selects outgroups from a single DIAMOND output file. """
    logging.debug(f"Processing DIAMOND file for OG {og_id}: {filepath.name}")
    try:
        # Define expected columns
        cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
        # Read TSV, handle potential bad lines
        df_hits = pd.read_csv(filepath, sep='\\t', header=None, names=cols,
                              low_memory=False, on_bad_lines='warn')

        # Convert numeric columns, coercing errors to NaN
        numeric_cols = ['evalue', 'bitscore', 'length', 'qlen', 'slen', 'pident',
                        'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send']
        for col in numeric_cols:
             if col in df_hits.columns:
                 df_hits[col] = pd.to_numeric(df_hits[col], errors='coerce')

        # Drop rows with NaN in essential numeric columns
        essential_cols = ['evalue', 'bitscore', 'length', 'qlen', 'slen']
        df_hits.dropna(subset=essential_cols, inplace=True)

        if df_hits.empty:
            logging.info(f"No valid hits found after initial cleaning for OG: {og_id}")
            return []

        # Classify hits based on keywords
        df_hits['classification'] = df_hits['sseqid'].apply(lambda x: get_outgroup_classification(x, tack_kw, eury_kw))

        # Calculate coverage (handle potential division by zero)
        df_hits['qcov'] = df_hits.apply(lambda row: row['length'] / row['qlen'] if row.get('qlen', 0) > 0 else 0, axis=1)
        df_hits['scov'] = df_hits.apply(lambda row: row['length'] / row['slen'] if row.get('slen', 0) > 0 else 0, axis=1)

        # Apply filters: E-value, Coverage, and Classification
        df_filtered = df_hits[
            (df_hits['evalue'] <= evalue_thr) &
            (df_hits['qcov'] >= min_cov) &
            (df_hits['scov'] >= min_cov) &
            (df_hits['classification'] != 'Unknown') # Filter out hits we couldn't classify
        ].copy()

        if df_filtered.empty:
            # Log why filtering failed (e.g., good hits existed but weren't classified)
            df_pre_tax_filter = df_hits[
                (df_hits['evalue'] <= evalue_thr) &
                (df_hits['qcov'] >= min_cov) &
                (df_hits['scov'] >= min_cov)
            ]
            if not df_pre_tax_filter.empty and (df_pre_tax_filter['classification'] == 'Unknown').all():
                 logging.info(f"Hits passed E-val/Cov for {og_id}, but none matched TACK/Eury keywords. No outgroups selected.")
            else:
                 logging.info(f"No hits passed E-val/Cov/Keyword filter for OG: {og_id}. No outgroups selected.")
            return [] # Return empty list if no hits pass filters

        # --- Selection Strategy ---
        # Sort by bitscore (descending) to prioritize best overall hits
        df_filtered = df_filtered.sort_values('bitscore', ascending=False)

        # Try to get the best TACK and best Eury hit first
        best_tack = df_filtered[df_filtered['classification'] == 'TACK'].head(1)
        best_eury = df_filtered[df_filtered['classification'] == 'Euryarchaeota'].head(1)

        # Combine best TACK/Eury, remove duplicates if the same hit was best for both (unlikely)
        potential_selection = pd.concat([best_tack, best_eury]).drop_duplicates(subset=['sseqid'])
        # Keep sorting by bitscore
        potential_selection = potential_selection.sort_values('bitscore', ascending=False)

        # Get the initial selection (up to max_outgroups)
        og_selection_ids = list(potential_selection['sseqid'])[:max_outgroups]

        # If we need more outgroups, fill with the next best overall hits
        num_needed = max_outgroups - len(og_selection_ids)
        if num_needed > 0:
            # Get remaining hits, excluding those already selected
            remaining_hits = df_filtered[~df_filtered['sseqid'].isin(og_selection_ids)]
            if not remaining_hits.empty:
                # Get unique IDs from the remaining hits, sorted by bitscore implicitly
                next_best_ids = list(remaining_hits['sseqid'].unique())[:num_needed]
                og_selection_ids.extend(next_best_ids)

        # Ensure final list has unique IDs and respects the maximum limit
        final_selection = list(dict.fromkeys(og_selection_ids))[:max_outgroups]
        logging.debug(f"Selected {len(final_selection)} outgroups for {og_id}: {final_selection}")
        return final_selection

    except pd.errors.EmptyDataError:
        logging.warning(f"Empty file or error reading columns for OG: {og_id} ({filepath.name}). Skipping.")
        return []
    except FileNotFoundError:
        logging.error(f"DIAMOND result file not found: {filepath}. Skipping.")
        return []
    except Exception as e:
        logging.exception(f"Processing file {filepath.name} for OG {og_id} failed: {e}. Skipping.")
        return []

def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Select outgroup sequences from DIAMOND search results based on keywords and scores.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_dir", required=True, type=Path,
                        help="Directory containing DIAMOND result files (TSV format).")
    parser.add_argument("-o", "--output_file", required=True, type=Path,
                        help="Path to save the output summary CSV file listing selected outgroups per OG.")
    parser.add_argument("--pattern", default="OG*_hits.tsv",
                        help="Glob pattern to find DIAMOND result files within the input directory.")
    parser.add_argument("--evalue", type=float, default=DEFAULT_EVALUE_THRESHOLD,
                        help="Maximum E-value threshold for filtering hits.")
    parser.add_argument("--min_cov", type=float, default=DEFAULT_MIN_COVERAGE,
                        help="Minimum query AND subject coverage threshold (0.0 to 1.0).")
    parser.add_argument("--max_outgroups", type=int, default=DEFAULT_MAX_OUTGROUPS_PER_OG,
                        help="Maximum number of outgroups to select per OG.")
    # Optional arguments for keyword lists (or load from file in future)
    # parser.add_argument("--tack_keywords_file", type=Path, help="File containing TACK keywords (one per line).")
    # parser.add_argument("--eury_keywords_file", type=Path, help="File containing Eury keywords (one per line).")
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Set the logging level.")

    return parser.parse_args()

# --- Main Execution Logic ---
def main():
    """ Main function to run the outgroup selection analysis. """
    args = parse_arguments()

    # Set logging level
    log_level_numeric = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.getLogger().setLevel(log_level_numeric)

    logging.info("--- Starting Outgroup Selection ---")
    logging.info(f"Input Directory: {args.input_dir}")
    logging.info(f"File Pattern: {args.pattern}")
    logging.info(f"Output File: {args.output_file}")
    logging.info(f"E-value Threshold: {args.evalue}")
    logging.info(f"Min Coverage: {args.min_cov}")
    logging.info(f"Max Outgroups per OG: {args.max_outgroups}")

    # --- Load Keywords (Using defaults for now) ---
    # Add logic here to load from files specified by args if implemented
    tack_keywords_set = set(DEFAULT_TACK_KEYWORDS)
    eury_keywords_set = set(DEFAULT_EURY_KEYWORDS)
    logging.info(f"Using {len(tack_keywords_set)} TACK keywords and {len(eury_keywords_set)} Eury keywords.")
    # ---

    # --- Validate Input Dir ---
    if not args.input_dir.is_dir():
        logging.error(f"Input directory not found: {args.input_dir}")
        sys.exit(1)

    # --- Find and Process Files ---
    all_selected_outgroups = defaultdict(list)
    search_path = args.input_dir / args.pattern
    # Use Path.glob for consistency
    diamond_files = sorted(list(args.input_dir.glob(args.pattern)))

    if not diamond_files:
        logging.error(f"No files found in '{args.input_dir}' matching pattern '{args.pattern}'")
        sys.exit(1)

    logging.info(f"Found {len(diamond_files)} potential DIAMOND result files...")
    processed_count = 0
    files_with_outgroups = 0

    for filepath in diamond_files:
        # Try to parse OG ID from filename (e.g., OG12345_...)
        og_match = re.match(r"(OG\\d+)", filepath.name)
        if og_match:
            og_id = og_match.group(1)
        else:
            # Fallback: Use filename stem if no OG pattern found
            og_id = filepath.stem
            logging.warning(f"Could not parse standard OG ID from '{filepath.name}', using '{og_id}' as identifier.")

        selected = process_diamond_file(
            filepath, og_id, tack_keywords_set, eury_keywords_set,
            args.evalue, args.min_cov, args.max_outgroups
        )
        all_selected_outgroups[og_id] = selected # Store even if empty list
        processed_count += 1
        if selected: # Only increment if list is not empty
             files_with_outgroups += 1

        # Log progress
        if processed_count % 100 == 0 or processed_count == len(diamond_files):
             logging.info(f"Processed {processed_count}/{len(diamond_files)} files...")
             # sys.stdout.flush() # Might not be needed with logging

    logging.info(f"\nFinished processing {processed_count} files.")
    logging.info(f"Found outgroups for {files_with_outgroups} OGs.")

    # --- Prepare and Save Output ---
    output_data = []
    # Sort OGs for consistent output order
    og_ids_sorted = sorted(all_selected_outgroups.keys())

    for og_id in og_ids_sorted:
        outgroups = all_selected_outgroups[og_id]
        row = {'OG_ID': og_id}
        # Create columns up to max_outgroups, fill with None if fewer selected
        for i in range(args.max_outgroups):
            col_name = f'Outgroup_{i+1}_sseqid'
            row[col_name] = outgroups[i] if i < len(outgroups) else None # Use None for missing
        output_data.append(row)

    if not output_data:
        logging.error("No outgroups selected for any OG. Output file will not be created.")
        return # Exit gracefully

    df_output = pd.DataFrame(output_data)
    # Ensure columns are in the correct order
    output_columns = ['OG_ID'] + [f'Outgroup_{i+1}_sseqid' for i in range(args.max_outgroups)]
    df_output = df_output.reindex(columns=output_columns)

    try:
        # Ensure output directory exists
        args.output_file.parent.mkdir(parents=True, exist_ok=True)
        # Save as CSV (can change separator if needed)
        df_output.to_csv(args.output_file, index=False, na_rep='NA') # Use NA for missing values
        logging.info(f"\nSuccessfully saved selected outgroups for {len(df_output)} OGs to: {args.output_file}")
    except Exception as e:
        logging.error(f"Could not save output file '{args.output_file}': {e}")

if __name__ == "__main__":
    main()
