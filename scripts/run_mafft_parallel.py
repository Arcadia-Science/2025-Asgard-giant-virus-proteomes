#!/usr/bin/env python3

"""
Runs MAFFT alignment in parallel on multiple input FASTA files using Python's
concurrent.futures module.

Skips files if the output alignment already exists and logs MAFFT's stderr
output for each job.
"""

import os
import sys
import glob
import subprocess
import logging
import time
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
# from multiprocessing import freeze_support # See note in main guard

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - [%(funcName)s] %(message)s',
                    stream=sys.stdout)

# --- Helper Function to Run MAFFT ---
def run_mafft_on_file(fasta_file_path: Path, output_dir: Path, log_dir: Path,
                      mafft_exe: str, mafft_args: list[str],
                      input_suffix: str, output_suffix: str) -> tuple[str, str]:
    """
    Runs MAFFT on a single input FASTA file.

    Args:
        fasta_file_path: Path to the input FASTA file.
        output_dir: Directory to save the output alignment.
        log_dir: Directory to save the MAFFT stderr log.
        mafft_exe: Path to the MAFFT executable.
        mafft_args: List of arguments for MAFFT (e.g., ["--auto", "--thread", "1"]).
        input_suffix: Expected suffix of input files (e.g., "_final.fasta").
        output_suffix: Suffix for output alignment files (e.g., "_final.aln").

    Returns:
        Tuple (input_filename, status_message).
        Status message is "Success", "Skipped (Output Exists)", or an error description.
    """
    base_name = fasta_file_path.name
    # Derive OG ID or base identifier by removing suffix
    if input_suffix and base_name.endswith(input_suffix):
        identifier = base_name[:-len(input_suffix)]
    else:
        identifier = fasta_file_path.stem # Fallback to stem if suffix doesn't match

    output_aln_path = output_dir / f"{identifier}{output_suffix}"
    output_log_path = log_dir / f"{identifier}_mafft.log"

    # Skip if output alignment already exists and is not empty
    if output_aln_path.exists() and output_aln_path.stat().st_size > 0:
        logging.debug(f"Output alignment {output_aln_path} already exists and is non-empty. Skipping.")
        return base_name, "Skipped (Output Exists)"
    elif output_aln_path.exists():
         logging.warning(f"Output alignment {output_aln_path} exists but is empty. Will overwrite.")

    # Construct MAFFT command
    command = [mafft_exe] + mafft_args + [str(fasta_file_path)]
    logging.debug(f"Running command for {identifier}: {' '.join(command)}")

    try:
        # Run MAFFT, redirect stdout to alignment file, stderr to log file
        with open(output_aln_path, 'w') as f_out, open(output_log_path, 'w') as f_err:
            # Set encoding if MAFFT output might contain non-ASCII characters
            process = subprocess.run(command, stdout=f_out, stderr=f_err, check=True, text=True, encoding='utf-8')

        # Check if output alignment is empty after running
        if not output_aln_path.stat().st_size > 0:
             error_info = "Unknown MAFFT issue (empty output)"
             try: # Try to read the log file for more info
                  with open(output_log_path, 'r', encoding='utf-8') as f_log_read:
                       log_content = f_log_read.read(500).strip() # Read first 500 chars
                       if log_content: error_info = f"MAFFT Error (see log): {log_content}..."
             except Exception: pass # Ignore errors reading the log file itself
             logging.warning(f"MAFFT created an empty alignment file for {base_name}. Check log: {output_log_path}")
             return base_name, error_info

        return base_name, "Success"

    except FileNotFoundError:
        # This error occurs if mafft_exe is not found
        logging.error(f"MAFFT executable '{mafft_exe}' not found.")
        # Return error status - this will likely happen for all jobs if path is wrong
        return base_name, f"Error: MAFFT not found ('{mafft_exe}')"
    except subprocess.CalledProcessError as e:
        # MAFFT exited with a non-zero status code
        logging.error(f"MAFFT failed for {base_name} with exit code {e.returncode}. Check log: {output_log_path}")
        return base_name, f"Error: MAFFT failed (Code {e.returncode})"
    except Exception as e:
        # Catch any other unexpected errors during subprocess execution
        logging.exception(f"An unexpected error occurred while running MAFFT for {base_name}: {e}")
        return base_name, f"Error: Unexpected ({type(e).__name__})"

def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run MAFFT alignment in parallel on multiple FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_dir", required=True, type=Path,
                        help="Input directory containing FASTA files.")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Output directory for aligned files.")
    parser.add_argument("-l", "--log_dir", required=True, type=Path,
                        help="Output directory for MAFFT stderr logs (.log).")
    parser.add_argument("--input_suffix", default=".fasta",
                        help="Suffix of input FASTA files to process (e.g., '.fasta', '_sequences.fa').")
    parser.add_argument("--output_suffix", default=".mafft.fa",
                        help="Suffix to use for output alignment files (default: .mafft.fa, for FastTree compatibility).")
    parser.add_argument("--mafft_exe", default="mafft",
                        help="Path to the MAFFT executable.")
    parser.add_argument("--mafft_args", default="--auto --thread 1",
                        help="Arguments to pass to MAFFT (as a single string, e.g., '--auto --quiet').")
    # Use os.cpu_count() to get a reasonable default, ensure at least 1 core
    default_cores = max(1, os.cpu_count() - 1 if os.cpu_count() and os.cpu_count() > 1 else 1)
    parser.add_argument("-n", "--num_cores", type=int, default=default_cores,
                        help="Number of parallel MAFFT jobs to run.")

    return parser.parse_args()

# --- Main Execution Logic ---
def main():
    """Main logic to run MAFFT alignments in parallel."""
    args = parse_arguments()
    start_time = time.time()

    # Validate MAFFT executable
    mafft_exe_path = shutil.which(args.mafft_exe)
    if not mafft_exe_path:
        logging.error(f"MAFFT executable '{args.mafft_exe}' not found in PATH or is not executable.")
        sys.exit(1)
    logging.info(f"Using MAFFT executable: {mafft_exe_path}")

    # Prepare MAFFT args list
    mafft_args_list = args.mafft_args.split()
    logging.info(f"Using MAFFT arguments: {mafft_args_list}")

    # Create output directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Output alignments will be saved in: {args.output_dir}")
    logging.info(f"Log files will be saved in: {args.log_dir}")

    # Find input FASTA files
    fasta_files = sorted(list(args.input_dir.glob(f"*{args.input_suffix}")))
    if not fasta_files:
        logging.error(f"No files matching '*{args.input_suffix}' found in {args.input_dir}. Stopping.")
        sys.exit(1)

    logging.info(f"Found {len(fasta_files)} FASTA files to align.")
    logging.info(f"Running up to {args.num_cores} MAFFT jobs in parallel...")

    success_count = 0
    skipped_count = 0
    error_count = 0
    results = []

    # Use ProcessPoolExecutor for parallel execution
    # max_workers should be <= number of physical cores for CPU-bound tasks like MAFFT
    with ProcessPoolExecutor(max_workers=args.num_cores) as executor:
        # Submit all jobs
        futures = {
            executor.submit(
                run_mafft_on_file,
                f, args.output_dir, args.log_dir,
                mafft_exe_path, mafft_args_list,
                args.input_suffix, args.output_suffix
            ): f
            for f in fasta_files
        }

        # Process completed jobs as they finish
        for i, future in enumerate(as_completed(futures)):
            input_filepath = futures[future]
            input_fname = input_filepath.name
            try:
                fname, status = future.result() # fname should match input_fname
                results.append({'file': fname, 'status': status})
                if status == "Success":
                    success_count += 1
                elif status == "Skipped (Output Exists)":
                    skipped_count += 1
                else: # Error status
                    error_count += 1
            except Exception as exc:
                logging.error(f"Job for {input_fname} generated an exception in the executor: {exc}")
                results.append({'file': input_fname, 'status': f"Error: Future Exception ({type(exc).__name__})"})
                error_count += 1

            # Log progress periodically
            processed_count = i + 1
            if processed_count % 50 == 0 or processed_count == len(fasta_files):
                 logging.info(f"Processed {processed_count}/{len(fasta_files)} files... (Success: {success_count}, Skipped: {skipped_count}, Errors: {error_count})")
                 # sys.stdout.flush() # Might not be needed with logging

    # --- Final Summary ---
    end_time = time.time()
    logging.info("--- MAFFT Alignment Summary ---")
    logging.info(f"Successfully aligned files: {success_count}")
    logging.info(f"Skipped (output existed): {skipped_count}")
    logging.info(f"Errors during alignment: {error_count}")
    logging.info(f"Total files processed/skipped: {success_count + skipped_count + error_count}/{len(fasta_files)}")
    logging.info(f"Aligned files saved in: {args.output_dir}")
    logging.info(f"Log files saved in: {args.log_dir}")
    logging.info(f"Total time: {end_time - start_time:.2f} seconds")
    logging.info("MAFFT alignment process completed.")

    if error_count > 0:
        logging.error(f"{error_count} jobs failed. Check logs in {args.log_dir} for details.")
        sys.exit(1) # Exit with error code if jobs failed

# --- Main execution Guard ---
if __name__ == '__main__':
    # freeze_support() # See note above regarding freeze_support
    import shutil # Needed for shutil.which check added in main()
    main()
