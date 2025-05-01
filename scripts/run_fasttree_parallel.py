#!/usr/bin/env python3

"""
Runs FastTree (preferably FastTreeMP) in parallel on multiple input alignment
files (FASTA format) using Python's concurrent.futures module.

Includes optional timeout per job and logs errors/timeouts. Skips files if the
output tree already exists and is non-empty.
"""

import os
import sys
import glob
import subprocess
import logging
import time
import argparse
import shutil # For checking executable path
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
# from multiprocessing import freeze_support # See note in main guard

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - [%(funcName)s] %(message)s',
                    stream=sys.stdout)

# --- Default Configuration ---
DEFAULT_FASTTREE_EXE = "FastTreeMP" # Prefer multi-threaded version
DEFAULT_N_CORES = max(1, os.cpu_count() - 1 if os.cpu_count() and os.cpu_count() > 1 else 1)
DEFAULT_FASTTREE_ARGS = ["-lg", "-gamma"] # Common args for protein, WAG/LG + CAT
DEFAULT_LOG_DIR_NAME = 'fasttree_logs'
DEFAULT_TIMEOUT_SECONDS = 3600 # 1 hour timeout, set to 0 or None to disable

# --- Function to Run FastTree ---
def run_fasttree_on_file(args_tuple: tuple) -> tuple[str, str]:
    """
    Runs FastTree on a single input trimmed alignment file with an optional timeout.
    Accepts a tuple of arguments for use with ProcessPoolExecutor.

    Args:
        args_tuple: Contains (trimmed_fasta_path, output_dir, log_dir,
                             fasttree_exe_path, fasttree_args_list, timeout_sec,
                             input_suffix, output_suffix)

    Returns:
        Tuple (input_filename, status_message).
    """
    # Unpack arguments
    trimmed_fasta_path, output_dir, log_dir, fasttree_exe_path, \
    fasttree_args_list, timeout_sec, input_suffix, output_suffix = args_tuple

    base_name = trimmed_fasta_path.name
    # Derive identifier by removing suffix
    if input_suffix and base_name.endswith(input_suffix):
        identifier = base_name[:-len(input_suffix)]
    else:
        identifier = trimmed_fasta_path.stem

    output_nwk_path = output_dir / f"{identifier}{output_suffix}"
    output_log_path = log_dir / f"{identifier}_fasttree.log"

    # Skip if output exists and is non-empty
    if output_nwk_path.exists() and output_nwk_path.stat().st_size > 0:
        logging.debug(f"Output tree {output_nwk_path} exists and is non-empty. Skipping.")
        return base_name, "Skipped (Output Exists)"
    elif output_nwk_path.exists():
        logging.warning(f"Output tree {output_nwk_path} exists but is empty. Will overwrite.")

    # Construct command
    command = [fasttree_exe_path] + fasttree_args_list
    logging.debug(f"Running command for {identifier}: {' '.join(command)} (Timeout: {timeout_sec}s)")

    try:
        # FastTree reads alignment from stdin, writes tree to stdout, logs to stderr
        with open(trimmed_fasta_path, 'r') as f_in, \
             open(output_nwk_path, 'w') as f_out, \
             open(output_log_path, 'w') as f_err:

            process = subprocess.run(command, stdin=f_in, stdout=f_out, stderr=f_err,
                                     check=True, text=True, encoding='utf-8', errors='ignore',
                                     timeout=timeout_sec) # Pass timeout value

        # Verify output file was created and is not empty
        if not output_nwk_path.exists() or not output_nwk_path.stat().st_size > 0:
             error_info = "FastTree Error (Output missing or empty)"
             try: # Try reading log for info
                  with open(output_log_path, 'r', encoding='utf-8') as f_log_read:
                       log_content = f_log_read.read(500).strip() # Read first 500 chars
                       if log_content: error_info = f"FastTree Info/Error (see log): {log_content}..."
             except Exception: pass
             # Attempt to remove empty file if it exists
             if output_nwk_path.exists():
                  try: output_nwk_path.unlink()
                  except OSError: pass
             logging.warning(f"FastTree ran for {base_name} but produced no/empty output. Check log: {output_log_path}")
             return base_name, error_info

        return base_name, "Success"

    except FileNotFoundError:
        return base_name, f"Error: FastTree not found at {fasttree_exe_path}"
    except subprocess.CalledProcessError as e:
        logging.error(f"FastTree failed for {base_name} with exit code {e.returncode}. Check log: {output_log_path}")
        if output_nwk_path.exists(): # Remove potentially incomplete output
             try: output_nwk_path.unlink()
             except OSError: pass
        return base_name, f"Error: FastTree failed (Code {e.returncode})"
    except subprocess.TimeoutExpired:
        logging.warning(f"FastTree timed out for {base_name} after {timeout_sec} seconds. Check log: {output_log_path}")
        if output_nwk_path.exists(): # Remove potentially incomplete output
             try: output_nwk_path.unlink()
             except OSError: pass
        return base_name, f"Error: Timeout ({timeout_sec}s)"
    except Exception as e:
        logging.exception(f"An unexpected error occurred while running FastTree for {base_name}: {e}")
        if output_nwk_path.exists(): # Remove potentially incomplete output
             try: output_nwk_path.unlink()
             except OSError: pass
        return base_name, f"Error: Unexpected ({type(e).__name__})"

def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run FastTree in parallel on multiple alignment files with timeout.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_dir", required=True, type=Path,
                        help="Directory containing input alignment files (FASTA format).")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Directory to save the output FastTree Newick files (.nwk).")
    parser.add_argument("-l", "--log_dir", default=None, type=Path,
                        help=f"Directory to save FastTree log files. Default: '<output_dir>/{DEFAULT_LOG_DIR_NAME}'")
    parser.add_argument("--input_suffix", default="_trimmed.fasta",
                        help="Suffix of input alignment files to process.")
    parser.add_argument("--output_suffix", default="_fasttree.nwk",
                        help="Suffix to use for output tree files.")
    parser.add_argument("-c", "--cores", type=int, default=DEFAULT_N_CORES,
                        help="Number of parallel FastTree jobs to run.")
    parser.add_argument("-f", "--fasttree_exe", default=DEFAULT_FASTTREE_EXE,
                        help="Path to the FastTree executable (FastTree or FastTreeMP).")
    parser.add_argument("--fasttree_args", default=" ".join(DEFAULT_FASTTREE_ARGS),
                        help="Arguments to pass to FastTree (as a single string, e.g., '-lg -gamma').")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_SECONDS,
                        help="Timeout in seconds for each FastTree job (0 or negative to disable).")
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Set the logging level.")

    return parser.parse_args()

# --- Main Execution Logic ---
def main():
    """ Main logic to run FastTree in parallel. """
    args = parse_arguments()
    start_time = time.time()

    # Set logging level
    log_level_numeric = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.getLogger().setLevel(log_level_numeric)

    logging.info("--- Starting FastTree Parallel Run ---")

    # Validate FastTree executable
    fasttree_exe_path = shutil.which(args.fasttree_exe)
    if not fasttree_exe_path:
        # Try plain 'FastTree' if default 'FastTreeMP' fails
        if args.fasttree_exe == DEFAULT_FASTTREE_EXE:
             logging.warning(f"Executable '{args.fasttree_exe}' not found, trying 'FastTree'...")
             fasttree_exe_path = shutil.which("FastTree")
        if not fasttree_exe_path:
             logging.error(f"FastTree executable ('{args.fasttree_exe}' or 'FastTree') not found in PATH.")
             sys.exit(1)
    logging.info(f"Using FastTree executable: {fasttree_exe_path}")

    # Prepare FastTree args list
    fasttree_args_list = args.fasttree_args.split()
    logging.info(f"Using FastTree arguments: {fasttree_args_list}")

    # Determine timeout value (None disables timeout in subprocess.run)
    timeout_value = args.timeout if args.timeout and args.timeout > 0 else None
    if timeout_value:
        logging.info(f"Setting timeout per job: {timeout_value} seconds")
    else:
         logging.info("No timeout set per job.")

    # Setup directories
    log_dir = args.log_dir if args.log_dir else args.output_dir / DEFAULT_LOG_DIR_NAME
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")
    logging.info(f"Log directory: {log_dir}")

    # Find input files
    fasta_files = sorted(list(args.input_dir.glob(f"*{args.input_suffix}")))
    if not fasta_files:
        logging.error(f"No files matching '*{args.input_suffix}' found in {args.input_dir}. Stopping.")
        sys.exit(1)
    logging.info(f"Found {len(fasta_files)} alignment files to process.")
    logging.info(f"Running up to {args.cores} FastTree jobs in parallel...")

    # --- Run in Parallel ---
    success_count = 0
    skipped_count = 0
    error_count = 0
    results_summary = []
    tasks = []
    for f in fasta_files:
        tasks.append((
            f, args.output_dir, log_dir, fasttree_exe_path,
            fasttree_args_list, timeout_value,
            args.input_suffix, args.output_suffix
        ))

    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = {executor.submit(run_fasttree_on_file, task): task[0] for task in tasks}
        for i, future in enumerate(as_completed(futures)):
            input_filepath = futures[future]
            input_fname = input_filepath.name
            processed_count = i + 1
            try:
                fname, status = future.result()
                results_summary.append({'file': fname, 'status': status})
                if status == "Success": success_count += 1
                elif status == "Skipped (Output Exists)": skipped_count += 1
                else:
                    error_count += 1
                    logging.warning(f"Job for {fname} failed with status: {status}") # Log errors immediately
            except Exception as exc:
                logging.error(f"Main loop error processing result for {input_fname}: {exc}")
                results_summary.append({'file': input_fname, 'status': f"Error: Future Exception ({type(exc).__name__})"})
                error_count += 1

            if processed_count % 50 == 0 or processed_count == len(fasta_files):
                 logging.info(f"Processed {processed_count}/{len(fasta_files)} files... (Success: {success_count}, Skipped: {skipped_count}, Errors: {error_count})")

    # --- Final Summary ---
    end_time = time.time()
    logging.info("--- FastTree Summary ---")
    logging.info(f"Successfully created trees: {success_count}")
    logging.info(f"Skipped (output existed): {skipped_count}")
    logging.info(f"Errors during tree building (incl. timeouts): {error_count}")
    logging.info(f"Total files processed/skipped: {success_count + skipped_count + error_count}/{len(fasta_files)}")
    logging.info(f"Tree files ({args.output_suffix}) saved in: {args.output_dir}")
    logging.info(f"Log files saved in: {log_dir}")

    if error_count > 0:
        logging.warning("--- Files with Errors/Timeouts ---")
        error_files_list = [r['file'] for r in results_summary if r['status'] not in ["Success", "Skipped (Output Exists)"]]
        for err_file in error_files_list[:20]: # Log first 20 errors
             logging.warning(f"  - {err_file}")
        if len(error_files_list) > 20: logging.warning("  ... (additional errors not listed)")
        logging.warning(f"Check corresponding logs in {log_dir} for details on the {error_count} error(s)/timeout(s).")

    logging.info(f"Total time: {end_time - start_time:.2f} seconds")
    logging.info("FastTree process completed.")

    if error_count > 0:
        sys.exit(1) # Exit with error code if jobs failed

# --- Main execution Guard ---
if __name__ == '__main__':
    # freeze_support() # See note in run_mafft_parallel.py
    main()
