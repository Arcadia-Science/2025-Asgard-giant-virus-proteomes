import os
from glob import glob
from joblib import Parallel, delayed
import argparse # ADDED: Import argparse

def run_fasttree(fasta_file, output_dir):
    """
    Runs FastTree on a single FASTA file.
    """
    try:
        base_name = os.path.basename(fasta_file).replace('.mafft.fa', '')
        output_file = os.path.join(output_dir, f"{base_name}.tree")
        
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        command = f"fasttree -nt -gtr < {fasta_file} > {output_file}"
        os.system(command)
        return f"Successfully processed {fasta_file}"
    except Exception as e:
        return f"Failed to process {fasta_file}: {e}"

def main():
    # ADDED: Set up argument parser
    parser = argparse.ArgumentParser(description="Run FastTree in parallel on aligned FASTA files.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing aligned FASTA files (e.g., from MAFFT).")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save the output .tree files.")
    parser.add_argument("-n", "--n_jobs", type=int, default=-1, help="Number of parallel jobs to run (-1 uses all available cores).")
    args = parser.parse_args()

    # MODIFIED: Use directories from command-line arguments
    aligned_files = glob(os.path.join(args.input_dir, '*.mafft.fa'))

    if not aligned_files:
        print(f"No aligned files found in {args.input_dir}")
        return

    print(f"Found {len(aligned_files)} aligned files to process.")
    print(f"Running on {args.n_jobs} cores.")

    # Run in parallel
    results = Parallel(n_jobs=args.n_jobs)(
        delayed(run_fasttree)(f, args.output_dir) for f in aligned_files
    )

    # Print results
    for r in results:
        print(r)

    print("FastTree parallel processing complete.")

if __name__ == "__main__":
    main()