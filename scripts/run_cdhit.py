#!/usr/bin/env python3

"""
Runs the CD-HIT clustering tool on an input FASTA file.

This script provides a wrapper around the cd-hit command-line tool,
allowing parameters like identity threshold, word size, memory, and threads
to be specified via arguments.
"""

import os
import sys
import argparse
import subprocess
import shlex
import logging
import shutil # To check if cd-hit executable exists

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def run_cdhit(input_fasta, output_fasta, identity_threshold, word_size, memory_limit_mb, threads, other_options=""):
    """
    Executes the cd-hit command using subprocess.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path for the output non-redundant FASTA file.
        identity_threshold (float): Sequence identity threshold (e.g., 0.9 for 90%).
        word_size (int): Word size (option -n in cd-hit).
        memory_limit_mb (int): Memory limit in MB (option -M in cd-hit, 0 for unlimited).
        threads (int): Number of threads (option -T in cd-hit).
        other_options (str): Any other cd-hit options as a single string (e.g., "-d 0").

    Returns:
        bool: True if cd-hit command executed successfully (exit code 0), False otherwise.
    """
    # Check if cd-hit executable is available
    cdhit_executable = shutil.which("cd-hit")
    if not cdhit_executable:
        logging.error("cd-hit executable not found in PATH. Please install CD-HIT.")
        return False

    # Construct the command
    # Use shlex.quote to handle paths with spaces safely
    cmd_parts = [
        cdhit_executable,
        "-i", shlex.quote(input_fasta),
        "-o", shlex.quote(output_fasta),
        "-c", str(identity_threshold),
        "-n", str(word_size),
        "-M", str(memory_limit_mb), # CD-HIT expects MB, 0 means unlimited
        "-T", str(threads)
    ]

    # Add other options if provided
    if other_options:
        # Use shlex.split to handle potentially complex options safely
        cmd_parts.extend(shlex.split(other_options))

    command_str = " ".join(cmd_parts) # For logging purposes
    logging.info(f"Running CD-HIT command:\n{command_str}")

    try:
        # Run the command, capture output and errors
        process = subprocess.run(command_str, shell=True, check=True, capture_output=True, text=True)
        logging.info("CD-HIT standard output:\n" + process.stdout)
        if process.stderr:
            logging.warning("CD-HIT standard error:\n" + process.stderr) # Log stderr as warning
        logging.info(f"CD-HIT finished successfully for {input_fasta}.")
        logging.info(f"Output non-redundant file: {output_fasta}")
        logging.info(f"Cluster information file: {output_fasta}.clstr")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"CD-HIT command failed with exit code {e.returncode}.")
        logging.error(f"Command: {e.cmd}")
        logging.error("Standard Output:\n" + e.stdout)
        logging.error("Standard Error:\n" + e.stderr)
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred while running CD-HIT: {e}")
        return False

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Run CD-HIT on a FASTA file.")
    parser.add_argument("-i", "--input_fasta", required=True,
                        help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output_fasta", required=True,
                        help="Path for the output non-redundant FASTA file.")
    parser.add_argument("-c", "--identity", type=float, default=0.9,
                        help="Sequence identity threshold (option -c, default: 0.9).")
    parser.add_argument("-n", "--word_size", type=int, choices=[2, 3, 4, 5], default=5,
                        help="Word size (option -n, default: 5). Must be 2, 3, 4, or 5.")
    parser.add_argument("-M", "--memory_mb", type=int, default=0,
                        help="Memory limit in MB (option -M, 0 for unlimited, default: 0).")
    parser.add_argument("-T", "--threads", type=int, default=1,
                        help="Number of threads (option -T, default: 1).")
    parser.add_argument("--other_options", type=str, default="-d 0",
                        help="String containing any other desired cd-hit options (e.g., '-d 0 -G 1') (default: '-d 0').")

    return parser.parse_args()

# --- Main Execution ---
if __name__ == "__main__":
    args = parse_arguments()

    # Validate input file existence
    if not os.path.isfile(args.input_fasta):
        logging.error(f"Input FASTA file not found: {args.input_fasta}")
        sys.exit(1)

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_fasta)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            logging.info(f"Created output directory: {output_dir}")
        except OSError as e:
            logging.error(f"Failed to create output directory {output_dir}: {e}")
            sys.exit(1)

    # Run CD-HIT
    success = run_cdhit(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        identity_threshold=args.identity,
        word_size=args.word_size,
        memory_limit_mb=args.memory_mb,
        threads=args.threads,
        other_options=args.other_options
    )

    if success:
        logging.info("Script finished successfully.")
        sys.exit(0)
    else:
        logging.error("Script finished with errors.")
        sys.exit(1)
