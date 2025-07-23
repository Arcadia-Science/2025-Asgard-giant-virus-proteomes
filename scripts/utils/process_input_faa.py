#!/usr/bin/env python3

"""
Processes protein FASTA files (.faa), typically from NCBI genome assemblies,
especially those with bracketed organism info common in viruses or MAGs.

1. Reads protein sequences from 'protein.faa' (or user-specified name) within
   each genome assembly subdirectory.
2. Parses the FASTA header to extract Protein ID and annotation. It specifically
   tries to extract and clean content within the first square brackets `[...]`
   to use as a 'Taxonomy/Source' field.
3. Cleans the main annotation string (part before brackets).
4. Determines annotation type (hypothetical/annotated).
5. Writes the sequences to a new FASTA file (one per genome) in an output
   directory, using a custom header format:
   >ProteinID|GenomeAssemblyID|BracketContent|AnnotationType|CleanedAnnotationName
"""

import os
import sys
import re
import argparse
import logging
import pandas as pd # For isnull check, can be removed if not strictly needed
from pathlib import Path
from Bio import SeqIO # Requires biopython

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# --- Helper Functions ---

def clean_annotation_name(annotation_string):
    """Cleans the protein annotation string (part before brackets)."""
    if pd.isnull(annotation_string) or not isinstance(annotation_string, str):
        return None
    # Take text before the first square bracket
    desc_part = annotation_string.split('[')[0].strip()
    if not desc_part:
        return None
    # Clean characters
    cleaned = re.sub(r'[\\s,;()\\[\\]{}:/\\\\]+', '_', desc_part)
    cleaned = re.sub(r'[^a-zA-Z0-9_.-]', '', cleaned) # Allow underscore, dot, hyphen
    cleaned = re.sub(r'_+', '_', cleaned).strip('_')
    if not cleaned: return None
    # Check against generic terms AFTER cleaning
    generic_terms = r'^(hypothetical_protein|unknown|predicted_protein|uncharacterized_protein|protein_of_unknown_function|possible_protein|orf|DUF.*)$'
    if re.match(generic_terms, cleaned, flags=re.IGNORECASE): return None
    return cleaned

def extract_clean_bracket_content(description):
    """Extracts and cleans content within the first square brackets `[...]`."""
    if not isinstance(description, str): return 'UnknownSource'
    match = re.search(r'\\[([^\\]]+)\\]', description) # Find text in first square brackets
    if match:
        raw_content = match.group(1).strip()
        # Clean characters (replace space, slash with underscore, allow dots/numbers/hyphen)
        content_clean = re.sub(r'[\\s/]+', '_', raw_content)
        content_clean = re.sub(r'[^a-zA-Z0-9_.-]', '', content_clean) # Allow underscore, dot, hyphen
        return content_clean if content_clean else 'UnknownSource' # Return cleaned or default
    return 'UnknownSource' # Default if no brackets found

def process_genome_faa(genome_dir_path: Path, output_dir_path: Path, input_faa_name: str, chars_per_line: int):
    """Processes a single genome's protein.faa file."""
    genome_id = genome_dir_path.name # Get GCA_... or other ID from directory name
    input_faa = genome_dir_path / input_faa_name
    output_fasta = output_dir_path / f"{genome_id}.fasta"

    if not input_faa.is_file():
        logging.warning(f"'{input_faa_name}' file not found in {genome_dir_path}. Skipping genome.")
        return 0, 0 # Return 0 processed, 0 written

    proteins_processed = 0
    proteins_written = 0
    records_to_write = []

    try:
        for record in SeqIO.parse(input_faa, "fasta"):
            proteins_processed += 1

            # Basic validation
            if not record.id or not record.seq:
                logging.warning(f"  Skipping invalid record in {input_faa.name} (ID: {record.id[:50]}...)")
                continue
            sequence = str(record.seq)
            if not sequence: # Check after conversion
                 logging.warning(f"  Skipping record with empty sequence in {input_faa.name} (ID: {record.id[:50]}...)")
                 continue

            # Parse header components
            protein_id = record.id
            description = record.description
            bracket_content = extract_clean_bracket_content(description)
            cleaned_name = clean_annotation_name(description) # Uses part before bracket

            # Determine annotation type based on description or cleaned name
            annotation_type = 'hypothetical' # Default
            generic_terms = r'(hypothetical|unknown|predicted|uncharacterized|domain_of_unknown_function)'
            if re.search(generic_terms, description, flags=re.IGNORECASE):
                annotation_type = 'hypothetical'
                cleaned_name = None # Ensure name is None if type is clearly hypothetical
            elif cleaned_name:
                annotation_type = 'annotated'
            else: # If cleaning failed but wasn't explicitly hypothetical
                annotation_type = 'hypothetical'
                cleaned_name = None

            # Create new header
            # Format: >ProteinID|GenomeAssemblyID|BracketContent|AnnotationType|CleanedAnnotationName
            name_field = cleaned_name if annotation_type == 'annotated' and cleaned_name else ''
            header_parts = [protein_id, genome_id, bracket_content, annotation_type, name_field]
            new_header = ">" + "|".join(str(p) for p in header_parts) # Ensure all parts are strings

            # Write record with new header and formatted sequence
            records_to_write.append(f"{new_header}\\n")
            for i in range(0, len(sequence), chars_per_line):
                records_to_write.append(sequence[i:i+chars_per_line] + '\\n')
            proteins_written += 1

        # Write all collected records to the output file
        if records_to_write:
            with open(output_fasta, 'w') as outfile:
                outfile.writelines(records_to_write)

        return proteins_processed, proteins_written

    except Exception as e:
        logging.error(f"Error processing {input_faa}: {e}")
        logging.exception("Traceback:") # Log traceback for unexpected errors
        return proteins_processed, 0 # Return processed count, but 0 written due to error

# --- Main Execution ---

def main():
    parser = argparse.ArgumentParser(description="Process NCBI protein FASTA files (.faa), parse headers/brackets, and standardize output FASTA.")
    parser.add_argument("-i", "--input_base_dir", required=True, type=Path,
                        help="Base directory containing genome assembly subdirectories.")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Directory to save the processed FASTA files.")
    parser.add_argument("--input_faa_name", default="protein.faa",
                        help="Name of the protein FASTA file within each genome subdirectory (default: protein.faa).")
    parser.add_argument("--line_length", type=int, default=60,
                        help="Number of characters per line in the output FASTA sequence (default: 60).")

    args = parser.parse_args()

    logging.info("Starting FASTA Processing and Header Standardization...")

    # --- Setup Directories ---
    input_base_path = Path(args.input_base_dir)
    output_path = Path(args.output_dir)
    if not input_base_path.is_dir():
        logging.error(f"Input base directory not found: {input_base_path}")
        sys.exit(1)
    try:
        output_path.mkdir(parents=True, exist_ok=True)
        logging.info(f"Output directory: {output_path}")
    except OSError as e:
        logging.error(f"Cannot create output directory {output_path}: {e}")
        sys.exit(1)

    # --- Process Each Genome Directory ---
    logging.info(f"Scanning for genome directories in: {input_base_path}")
    # Assume genome directories start with GCA_ or GCF_ or are just directories
    genome_dirs = [d for d in input_base_path.iterdir() if d.is_dir()]
    logging.info(f"Found {len(genome_dirs)} potential genome directories.")

    total_proteins_processed = 0
    total_proteins_written = 0
    genomes_processed_count = 0
    genomes_skipped_count = 0
    genomes_with_errors = 0

    for genome_dir in sorted(genome_dirs):
        logging.info(f"\\nProcessing Genome Directory: {genome_dir.name}")
        processed, written = process_genome_faa(genome_dir, output_path, args.input_faa_name, args.line_length)

        if processed == 0 and written == 0: # Indicates file not found
            genomes_skipped_count += 1
        elif written == 0 and processed > 0: # Indicates error during processing
             genomes_with_errors += 1
             genomes_processed_count += 1 # Still count as processed if attempted
             total_proteins_processed += processed # Add processed count even on error
        else: # Success
            genomes_processed_count += 1
            total_proteins_processed += processed
            total_proteins_written += written
            logging.info(f"  Finished. Processed {processed:,} proteins. Wrote {written:,} to: {genome_dir.name}.fasta")

    # --- Final Summary ---
    logging.info("\\n" + "="*60)
    logging.info("Processing Complete.")
    logging.info(f"Genome directories scanned: {len(genome_dirs)}")
    logging.info(f"Genomes processed (attempted): {genomes_processed_count}")
    logging.info(f"Genomes skipped ('{args.input_faa_name}' not found): {genomes_skipped_count}")
    logging.info(f"Genomes with processing errors: {genomes_with_errors}")
    logging.info(f"Total proteins read across all files: {total_proteins_processed:,}")
    logging.info(f"Total proteins written: {total_proteins_written:,}")
    logging.info(f"Output FASTA files are in: {args.output_dir}")
    logging.info("="*60)

    if genomes_with_errors > 0:
         logging.warning("Processing finished with errors for some genomes. Check logs.")
         sys.exit(1)

if __name__ == "__main__":
    main()
