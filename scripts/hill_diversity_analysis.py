"""
Calculates Hill diversity profiles for orthogroups based on feature abundances.

This script takes a CSV file where rows represent orthogroups and columns represent
the abundance of different features (e.g., domain counts, species counts) within
each orthogroup. It calculates the Hill diversity for a range of `q` values (from 0 to 3)
for each orthogroup.

Hill diversity is a unified framework for measuring biodiversity. The parameter `q`
determines the sensitivity of the metric to the rarity of features:
- q=0: Richness (the number of unique features).
- q=1: The exponential of Shannon entropy (weights features by their frequency).
- q=2: The inverse Simpson index (gives more weight to common features).

The output is a long-format CSV file containing the orthogroup, the `q` value, and the
calculated diversity, suitable for plotting diversity profiles.

Example:
    python scripts/hill_diversity_analysis.py \\
        -i data/analysis/orthogroup_domain_counts.csv \\
        -o results/diversity_profiles/domain_hill_diversity.csv
"""
import pandas as pd
import numpy as np
import argparse
import logging

# --- Setup Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def hill_diversity(p, q):
    """
    Calculates Hill diversity for a given probability vector p and order q.
    """
    # Filter out zero probabilities to avoid issues with log(0)
    p = p[p > 0]
    if q == 1:
        # Special case for q=1, which is the exponential of Shannon entropy
        return np.exp(-np.sum(p * np.log(p)))
    else:
        return np.sum(p**q) ** (1 / (1 - q))


def main():
    """
    Main function to calculate Hill diversity metrics from orthogroup data
    and save the results to a CSV file.
    """
    parser = argparse.ArgumentParser(
        description="Calculate Hill diversity from orthogroup domain abundance data."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input CSV file. Expects Orthogroups as rows and domain/species counts as columns.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to save the output CSV file containing calculated diversity metrics.",
    )
    args = parser.parse_args()

    logging.info(f"Reading input data from: {args.input}")
    try:
        df = pd.read_csv(args.input)
        # Assume the first column is the orthogroup identifier
        df = df.set_index(df.columns[0])
    except FileNotFoundError:
        logging.error(f"Input file not found: {args.input}")
        return
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        return

    # Define the range of q values for the diversity profile
    q_values = np.linspace(0, 3, 50)
    results = []

    logging.info(f"Calculating Hill diversity for {len(df)} orthogroups...")
    for orthogroup, row in df.iterrows():
        # Convert counts to proportions
        total_count = np.sum(row.values)
        if total_count == 0:
            continue  # Skip orthogroups with no domains/species

        proportions = row.values / total_count

        for q in q_values:
            diversity = hill_diversity(proportions, q)
            results.append({"Orthogroup": orthogroup, "q": q, "Diversity": diversity})

    if not results:
        logging.warning("No diversity results were calculated. Check input data.")
        return

    results_df = pd.DataFrame(results)

    logging.info(f"Saving diversity metrics to: {args.output}")
    try:
        results_df.to_csv(args.output, index=False, float_format="%.5f")
        logging.info("Analysis complete.")
    except IOError as e:
        logging.error(f"Failed to write output file: {e}")


if __name__ == "__main__":
    main()
