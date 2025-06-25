import os
import sys
import glob
import subprocess
import logging
import time
import argparse
import shutil
from pathlib import Path
from joblib import Parallel, delayed

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - [%(funcName)s] %(message)s',
                    stream=sys.stdout)

# --- Helper Function to Run FastTree ---
def run_fasttree_on_file(aligned_fasta_path: Path,
                         output_dir: Path,
                         log_dir: Path,
                         fasttree_exe: str,
                         fasttree_args: list[str],
                         input_suffix: str,
                         output_suffix: str) -> tuple[str, str]:
    """
    Runs FastTree on a single input aligned FASTA file.

    Args:
        aligned_fasta_path: Path to the input aligned FASTA file.
        output_dir: Directory to save the output tree file.
        log_dir: Directory to save the FastTree stderr log.
        fasttree_exe: Path to the FastTree executable.
        fasttree_args: List of arguments for FastTree (e.g., ["-nt", "-gtr"]).
        input_suffix: Expected suffix of input files (e.g., ".mafft.fa", ".aln").
        output_suffix: Suffix for output tree files (e.g., ".tree").

    Returns:
        Tuple (input_filename, status_message).
        Status message is "Success", "Skipped (Output Exists)", or an error description.
    """
    base_name = aligned_fasta_path.name
    # Derive identifier by removing input_suffix
    if input_suffix and base_name.endswith(input_suffix):
        identifier = base_name[:-len(input_suffix)]
    else:
        identifier = aligned_fasta_path.stem # Fallback to stem

    output_tree_path = output_dir / f"{identifier}{output_suffix}"
    output_log_path = log_dir / f"{identifier}_fasttree.log"

    # Skip if output tree already exists and is not empty
    if output_tree_path.exists() and output_tree_path.stat().st_size > 0:
        logging.debug(f"Output tree {output_tree_path} already exists and is non-empty. Skipping.")
        return base_name, "Skipped (Output Exists)"
    elif output_tree_path.exists():
        logging.warning(f"Output tree {output_tree_path} exists but is empty. Will overwrite.")

    # Construct FastTree command: fasttree [args] < input > output
    # Stdout is redirected to the tree file, Stderr to the log file.
    command = [fasttree_exe] + fasttree_args
    logging.debug(f"Running command for {identifier}: {' '.join(command)} < {aligned_fasta_path} > {output_tree_path}")

    try:
        with open(aligned_fasta_path, 'r') as stdin_file, \
             open(output_tree_path, 'w') as stdout_file, \
             open(output_log_path, 'w') as stderr_file:

            process = subprocess.run(command, stdin=stdin_file, stdout=stdout_file, stderr=stderr_file,
                                     check=True, text=True, encoding='utf-8')

        # Check if output tree is empty after running
        if not output_tree_path.stat().st_size > 0:
            error_info = "Unknown FastTree issue (empty output tree)"
            try: # Try to read the log file for more info
                with open(output_log_path, 'r', encoding='utf-8') as f_log_read:
                    log_content = f_log_read.read(500).strip() # Read first 500 chars
                    if log_content: error_info = f"FastTree Error (see log): {log_content}..."
            except Exception: pass
            logging.warning(f"FastTree created an empty tree file for {base_name}. Check log: {output_log_path}")
            return base_name, error_info
        
        return base_name, "Success"

    except FileNotFoundError:
        logging.error(f"FastTree executable '{fasttree_exe}' not found.")
        return base_name, f"Error: FastTree not found ('{fasttree_exe}')"
    except subprocess.CalledProcessError as e:
        logging.error(f"FastTree failed for {base_name} with exit code {e.returncode}. Check log: {output_log_path}")
        return base_name, f"Error: FastTree failed (Code {e.returncode})"
    except Exception as e:
        logging.exception(f"An unexpected error occurred while running FastTree for {base_name}: {e}")
        return base_name, f"Error: Unexpected ({type(e).__name__})"

def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run FastTree in parallel on multiple aligned FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input_dir", required=True, type=Path,
                        help="Directory containing aligned FASTA files.")
    parser.add_argument("-o", "--output_dir", required=True, type=Path,
                        help="Directory to save the output tree files.")
    parser.add_argument("-l", "--log_dir", required=True, type=Path,
                        help="Directory to save FastTree stderr logs.")
    parser.add_argument("--input_suffix", default=".mafft.fa",
                        help="Suffix of input aligned FASTA files to process (e.g., '.aln', '.mafft.fa').")
    parser.add_argument("--output_suffix", default=".tree",
                        help="Suffix for output tree files.")
    parser.add_argument("--fasttree_exe", default="fasttree",
                        help="Path to the FastTree executable.")
    parser.add_argument("--fasttree_args", default="-nt -gtr",
                        help="Arguments to pass to FastTree (as a single string, e.g., '-nt -lg').")
    default_cores = max(1, os.cpu_count() - 2 if os.cpu_count() and os.cpu_count() > 2 else 1)
    parser.add_argument("-n", "--n_jobs", type=int, default=default_cores,
                        help="Number of parallel FastTree jobs to run.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    start_time = time.time()

    fasttree_exe_path = shutil.which(args.fasttree_exe)
    if not fasttree_exe_path:
        logging.error(f"FastTree executable '{args.fasttree_exe}' not found in PATH or is not executable.")
        sys.exit(1)
    logging.info(f"Using FastTree executable: {fasttree_exe_path}")

    fasttree_args_list = args.fasttree_args.split()
    logging.info(f"Using FastTree arguments: {fasttree_args_list}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Output trees will be saved in: {args.output_dir}")
    logging.info(f"Log files will be saved in: {args.log_dir}")

    aligned_files = sorted(list(args.input_dir.glob(f"*{args.input_suffix}")))
    if not aligned_files:
        logging.error(f"No files matching '*{args.input_suffix}' found in {args.input_dir}. Stopping.")
        sys.exit(1)

    logging.info(f"Found {len(aligned_files)} aligned files to process for tree building.")
    logging.info(f"Running up to {args.n_jobs} FastTree jobs in parallel...")

    success_count = 0
    skipped_count = 0
    error_count = 0

    results_log = []

    with Parallel(n_jobs=args.n_jobs) as parallel:
        futures = parallel(
            delayed(run_fasttree_on_file)(
                f, args.output_dir, args.log_dir,
                fasttree_exe_path, fasttree_args_list,
                args.input_suffix, args.output_suffix
            ) for f in aligned_files
        )

        for i, (fname, status) in enumerate(futures):
            results_log.append({'file': fname, 'status': status})
            if status == "Success":
                success_count += 1
            elif status == "Skipped (Output Exists)":
                skipped_count += 1
            else: # Error status
                error_count += 1

            processed_count = i + 1
            if processed_count % 50 == 0 or processed_count == len(aligned_files):
                 logging.info(f"Processed {processed_count}/{len(aligned_files)} files... (Success: {success_count}, Skipped: {skipped_count}, Errors: {error_count})")

    end_time = time.time()
    logging.info("--- FastTree Processing Summary ---")
    logging.info(f"Successfully built trees for: {success_count} files")
    logging.info(f"Skipped (output existed): {skipped_count}")
    logging.info(f"Errors during tree building: {error_count}")
    logging.info(f"Total files processed/skipped: {success_count + skipped_count + error_count}/{len(aligned_files)}")
    logging.info(f"Tree files saved in: {args.output_dir}")
    logging.info(f"Log files saved in: {args.log_dir}")
    logging.info(f"Total time: {end_time - start_time:.2f} seconds")
    logging.info("FastTree parallel processing complete.")

    if error_count > 0:
        logging.error(f"{error_count} jobs failed. Check logs in {args.log_dir} for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()