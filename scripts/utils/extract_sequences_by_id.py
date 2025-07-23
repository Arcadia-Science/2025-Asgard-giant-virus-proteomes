#!/usr/bin/env python3

"""
Extracts sequences from a large reference FASTA file based on sequence IDs
found in multiple 'hit' files (typically TSV/CSV format from tools like DIAMOND).

Assumes each hit file corresponds to an Orthologous Group (OG) or query,
and contains a column with the sequence IDs (sseqid) to extract from the
reference FASTA. Creates one output FASTA file per input hit file.

Includes logic to handle potential discrepancies between the ID format in the
hit file and the full header format in the reference FASTA file index.
"""

import os
import sys
import glob
import pandas as pd
from Bio import SeqIO # Requires biopython
from Bio.SeqRecord import SeqRecord
import argparse
import logging
from pathlib import Path

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] [%(funcName)s] %(message)s', stream=sys.stdout)

# --- Helper Functions ---

def index_fasta(fasta_file: Path) -> SeqIO.index | None:
    """ Indexes a FASTA file for fast sequence retrieval using SeqIO.index. """
    logging.info(f"Indexing reference FASTA file: {fasta_file} ... (may take a moment)")
    try:
        # 'r' mode is implicit, specify alphabet if needed, e.g., alphabet=Alphabet.generic_protein
        fasta_index = SeqIO.index(str(fasta_file), "fasta")
        logging.info(f"Successfully indexed {len(fasta_index):,} sequences.")
        return fasta_index
    except FileNotFoundError:
        logging.error(f"Reference FASTA file not found: {fasta_file}")
        return None
    except Exception as e:
        logging.error(f"Error indexing FASTA file {fasta_file}: {e}")
        return None

def build_lookup_map(fasta_index: SeqIO.index) -> dict[str, str]:
    """
    Builds a map from a simplified key (e.g., part before first '|')
    to the actual full key used by the Biopython FASTA index.

    This helps find sequences even if the hit file ID doesn't exactly match
    the full FASTA header indexed by Biopython.
    """
    logging.info("Building lookup key map from FASTA index...")
    lookup_key_map = {}
    duplicates = 0
    map_keys_logged = 0

    total_keys = len(fasta_index)
    logging.info(f"Processing {total_keys:,} keys from FASTA index...")
    processed_keys = 0

    for actual_biopython_key in fasta_index.keys():
        processed_keys += 1
        if processed_keys % 100000 == 0:
             logging.info(f"  Mapped {processed_keys:,}/{total_keys:,} keys...")

        # --- Key Derivation Logic (ADJUST IF NEEDED) ---
        # Default: Use part before the first pipe '|' as the lookup key.
        # Modify this if your hit file IDs match a different part of the header.
        lookup_key = str(actual_biopython_key).split('|', 1)[0].strip()
        # --- End Key Derivation Logic ---

        if not lookup_key: # Skip if derivation results in empty key
            logging.debug(f"Skipping empty lookup key derived from actual key: '{actual_biopython_key}'")
            continue

        if lookup_key in lookup_key_map:
            # Handle cases where multiple FASTA entries might map to the same lookup key.
            # Keep the first one encountered.
            if duplicates < 5: # Log only the first few duplicates
                 logging.warning(f"Lookup key '{lookup_key}' derived from multiple FASTA index keys. Keeping mapping to first encountered: '{lookup_key_map[lookup_key]}'. Ignoring mapping from '{actual_biopython_key}'.")
            elif duplicates == 5:
                 logging.warning("... (suppressing further duplicate key warnings)")
            duplicates += 1
        else:
            lookup_key_map[lookup_key] = actual_biopython_key
            # Log first few mappings for debugging
            if map_keys_logged < 5:
                 logging.debug(f"Map: lookup_key='{lookup_key}' -> actual_key='{actual_biopython_key}'")
                 map_keys_logged +=1
            elif map_keys_logged == 5:
                 logging.debug("... (suppressing further map key logging)")
                 map_keys_logged += 1

    logging.info(f"Built map with {len(lookup_key_map):,} unique lookup keys. Encountered {duplicates:,} potential duplicate mappings.")
    return lookup_key_map


def extract_sequences_for_og(hits_file: Path, fasta_index: SeqIO.index, lookup_key_map: dict,
                             output_dir: Path, hits_column: str, header_column: str | None):
    """
    Extracts sequences listed in a hits file using the lookup map.

    Args:
        hits_file: Path to the input hits file (TSV/CSV).
        fasta_index: Indexed reference FASTA file object.
        lookup_key_map: Dictionary mapping simplified IDs to full FASTA index keys.
        output_dir: Directory to save the output FASTA file.
        hits_column: Name of the column in hits_file containing sequence IDs.
        header_column: Optional name of column containing full headers to use in output.
    """
    try:
        # Try to infer OG ID from filename (e.g., OG000123_best_euk_hits.txt -> OG000123)
        # Customize this regex if filename patterns differ
        match = re.match(r"^(OG\\d+)", hits_file.name)
        if match:
            og_id = match.group(1)
        else:
            og_id = hits_file.stem # Fallback to filename without extension
        output_fasta = output_dir / f"{og_id}_extracted_sequences.fasta"
        logging.info(f"Extracting sequences for {og_id} based on {hits_file.name}...")

        extracted_records = []
        not_found_count = 0
        processed_count = 0
        ids_in_hits_file = set()

        # Read hits file - try tab first, then comma
        try:
            hits_df = pd.read_csv(hits_file, sep='\\t')
            if hits_column not in hits_df.columns:
                logging.debug(f"'{hits_column}' not found with tab separator, trying comma for {hits_file.name}...")
                hits_df = pd.read_csv(hits_file, sep=',')
        except Exception as read_e:
             logging.error(f"Failed to read hits file {hits_file.name} with tab or comma separator: {read_e}")
             return # Skip this file

        # Validate required columns
        required_cols = [hits_column]
        if header_column: required_cols.append(header_column)
        if not all(col in hits_df.columns for col in required_cols):
             logging.error(f"Hits file {hits_file.name} missing required column(s): {', '.join(required_cols)}. Found: {list(hits_df.columns)}")
             return # Skip this file

        if hits_df.empty:
            logging.warning(f"Hits file {hits_file.name} is empty.")
            # Create empty output file to mark as processed
            output_fasta.touch()
            return

        # Process hits
        for index, row in hits_df.iterrows():
            processed_count += 1
            sseqid_raw = row[hits_column]
            if pd.isna(sseqid_raw):
                 logging.debug(f"Skipping row {index+2} in {hits_file.name} due to missing '{hits_column}'.")
                 continue
            sseqid = str(sseqid_raw).strip()
            ids_in_hits_file.add(sseqid) # Track unique IDs requested

            # --- Key Lookup Logic ---
            # 1. Try direct lookup using the sseqid as provided
            actual_biopython_key = lookup_key_map.get(sseqid)
            # 2. If not found, derive simplified key (e.g., before '|') and try again
            if not actual_biopython_key:
                lookup_key_derived = sseqid.split('|', 1)[0].strip()
                if lookup_key_derived != sseqid: # Only try derived key if different
                     actual_biopython_key = lookup_key_map.get(lookup_key_derived)
                     if actual_biopython_key:
                          logging.debug(f"Found match for '{sseqid}' using derived key '{lookup_key_derived}' -> '{actual_biopython_key}'")
            # --- End Key Lookup Logic ---

            if actual_biopython_key:
                try:
                    # Retrieve the record using the key found in the map
                    original_record = fasta_index[actual_biopython_key]

                    # Determine header for output FASTA
                    output_header_id = sseqid # Default to the ID from the hits file
                    output_description = ""
                    if header_column and header_column in row and pd.notna(row[header_column]):
                         # Use the full header from the specified column if available
                         full_header_from_file = str(row[header_column]).strip()
                         # Ensure the ID part matches sseqid, otherwise just use sseqid
                         if full_header_from_file.startswith(sseqid):
                              output_header_id = full_header_from_file
                         else:
                              output_description = full_header_from_file # Put mismatching header in description

                    # Create a new SeqRecord
                    new_record = SeqRecord(
                        seq=original_record.seq,
                        id=output_header_id,
                        description=output_description
                    )
                    extracted_records.append(new_record)

                except KeyError:
                    logging.warning(f"Mapped key '{actual_biopython_key}' for ID '{sseqid}' not found in FASTA index unexpectedly! Skipping.")
                    not_found_count += 1
                except Exception as e:
                     logging.error(f"Error processing sequence '{sseqid}' after finding map key '{actual_biopython_key}': {e}")
                     not_found_count += 1
            else:
                # ID from hits file (and its derived key) was not found in our map
                logging.warning(f"ID '{sseqid}' (from hits file {hits_file.name}) could not be mapped to any key in the reference FASTA index. Skipping.")
                not_found_count += 1

        # Write the extracted records to the output FASTA file
        if extracted_records:
            try:
                # Ensure output directory exists
                output_dir.mkdir(parents=True, exist_ok=True)
                with open(output_fasta, 'w') as outfile:
                    SeqIO.write(extracted_records, outfile, "fasta")
                logging.info(f"Saved {len(extracted_records)} sequences to {output_fasta}. ({not_found_count} sequences not found/skipped from {len(ids_in_hits_file)} unique IDs in hits file).")
            except Exception as e:
                logging.error(f"Error writing output FASTA file {output_fasta}: {e}")
        else:
            logging.warning(f"No sequences were successfully extracted for {og_id} based on {hits_file.name}. Output file may be empty or not created.")
            # Create empty file to mark as processed
            output_fasta.touch()

    except pd.errors.EmptyDataError:
        logging.warning(f"Hits file {hits_file.name} is empty or could not be parsed. Skipping.")
        (output_dir / f"{og_id}_extracted_sequences.fasta").touch() # Create empty file
    except FileNotFoundError:
        logging.error(f"Hits file not found: {hits_file}")
    except Exception as e:
        logging.exception(f"Failed to process hits file {hits_file.name}: {e}")


def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Extract sequences from a reference FASTA based on IDs listed in multiple hit files.")
    parser.add_argument("-r", "--ref_fasta", required=True, type=Path,
                        help="Path to the large reference FASTA file containing all sequences.")
    parser.add_argument("-d", "--hits_dir", required=True, type=Path,
                        help="Directory containing the hit files (e.g., *_best_euk_hits.txt).")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Directory to save the per-OG extracted FASTA files.")
    parser.add_argument("--hits_suffix", default="_best_euk_hits.txt",
                        help="Suffix of the hit files to process (default: _best_euk_hits.txt).")
    parser.add_argument("--hits_column", default="sseqid",
                        help="Name of the column in hit files containing sequence IDs to extract (default: sseqid).")
    parser.add_argument("--header_column", default="full_header",
                        help="Optional: Name of column in hit files containing full headers to use in output FASTA (default: full_header). If None or column missing, uses the ID from --hits_column.")

    return parser.parse_args()

# --- Main Execution ---
if __name__ == "__main__":
    args = parse_arguments()

    # --- Validate paths ---
    if not args.ref_fasta.is_file():
        logging.error(f"Input reference FASTA not found: {args.ref_fasta}")
        sys.exit(1)
    if not args.hits_dir.is_dir():
        logging.error(f"Input hits directory not found: {args.hits_dir}")
        sys.exit(1)
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logging.error(f"Cannot create output directory '{args.output_dir}': {e}")
        sys.exit(1)

    # --- Index FASTA and Build Map ---
    fasta_index = index_fasta(args.ref_fasta)
    if fasta_index is None:
        logging.error("Failed to index reference FASTA file. Exiting.")
        sys.exit(1)
    lookup_key_map = build_lookup_map(fasta_index)
    if not lookup_key_map:
        logging.error("Failed to build lookup key map from FASTA index. Exiting.")
        if hasattr(fasta_index, 'close'): fasta_index.close() # Close index if possible
        sys.exit(1)

    # --- Find and Process Hit Files ---
    hits_files = sorted(list(args.hits_dir.glob(f"*{args.hits_suffix}")))
    if not hits_files:
        logging.error(f"No hit files matching '*{args.hits_suffix}' found in {args.hits_dir}. Exiting.")
        if hasattr(fasta_index, 'close'): fasta_index.close()
        sys.exit(1)

    logging.info(f"Found {len(hits_files)} hit files to process.")

    # --- Process each hits file ---
    files_processed = 0
    for hits_file in hits_files:
        files_processed += 1
        extract_sequences_for_og(
            hits_file,
            fasta_index,
            lookup_key_map,
            args.output_dir,
            args.hits_column,
            args.header_column
        )
        # Optional: Flush stdout buffer if running interactively and want immediate feedback
        # sys.stdout.flush()

    # --- Cleanup ---
    if hasattr(fasta_index, 'close'):
         fasta_index.close() # Close the index file handle if SeqIO.index opened one
    logging.info(f"Finished processing {files_processed} hit files.")
    logging.info(f"Output FASTA files are in: {args.output_dir}")
