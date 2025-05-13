#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script to calculate Hill's diversity (q=1) for phylogenetic trees and multiple sequence alignments (MSAs)
# Based on the logic from the provided R script hills_diversity_tree_msa.R
# Combined into a single script for easier use.

import math
import numpy as np
import os
import glob
from collections import Counter
from Bio import Phylo
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
import pandas as pd

# Ensure necessary libraries are installed:
# pip install biopython numpy pandas

def hill_diversity_tree_py(tree: Tree, epsilon: float = 1e-10) -> dict:
    """
    Calculates Hill's diversity (q=1, exponential of Shannon entropy) for a phylogenetic tree.
    Also provides normalized measures.

    Args:
        tree: A Biopython Phylo.BaseTree.Tree object.
        epsilon: A small scalar added to branch lengths to avoid log(0) issues.

    Returns:
        A dictionary containing:
        - Hill_diversity: Hill's diversity of order 1 (exp(Shannon entropy of branch lengths)).
        - Norm_by_tips: Hill_diversity normalized by the number of tips.
        - Norm_by_PD: Hill_diversity normalized by Faith's Phylogenetic Diversity.
        - n_tips: The number of tips in the tree.
        - entropy: The Shannon entropy of branch length proportions.
    """
    branches = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]

    if not branches:
        # print("Warning: No branch lengths found in the tree. Cannot calculate diversity.")
        return {
            "Hill_diversity": None,
            "Norm_by_tips": None,
            "Norm_by_PD": None,
            "n_tips": len(tree.get_terminals()),
            "entropy": None,
        }

    branches_plus_epsilon = [b + epsilon for b in branches]

    total_length_plus_epsilon = sum(branches_plus_epsilon)
    if total_length_plus_epsilon == 0:
         print("Warning: Sum of branch lengths is zero even with epsilon. Cannot calculate proportions.")
         return {
            "Hill_diversity": None,
            "Norm_by_tips": None,
            "Norm_by_PD": None,
            "n_tips": len(tree.get_terminals()),
            "entropy": None,
        }

    p_branch = [b / total_length_plus_epsilon for b in branches_plus_epsilon]

    entropy = -sum(p * math.log(p) for p in p_branch if p > 0)

    hill_div = math.exp(entropy)

    n_tips = len(tree.get_terminals())
    original_branches_sum = sum([clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None])

    norm_by_tips = hill_div / n_tips if n_tips > 0 else None
    norm_by_pd = hill_div / original_branches_sum if original_branches_sum is not None and original_branches_sum > 0 else None

    return {
        "Hill_diversity": hill_div,
        "Norm_by_tips": norm_by_tips,
        "Norm_by_PD": norm_by_pd,
        "n_tips": n_tips,
        "entropy": entropy,
    }

def hill_diversity_msa_py(msa: MultipleSeqAlignment) -> dict:
    """
    Calculates Hill's diversity (q=1, exponential of Shannon entropy) for each
    column of a multiple sequence alignment (MSA), and provides mean and normalized mean.

    Args:
        msa: A Biopython Bio.Align.MultipleSeqAlignment object.

    Returns:
        A dictionary containing:
        - column_hill: A list of Hill's diversity values for each column.
        - mean_diversity: The mean Hill's diversity across all columns.
        - normalized_diversity: The mean diversity normalized by the number of possible amino acids (20).
    """
    if not msa or len(msa) == 0:
        # print("Warning: Input MSA is empty. Cannot calculate diversity.")
        return {
            "column_hill": [],
            "mean_diversity": None,
            "normalized_diversity": None
        }

    S = 20.0 # Number of possible states (amino acids)

    column_hill = []
    # Correct way to iterate through columns of a MultipleSeqAlignment
    alignment_length = msa.get_alignment_length()
    if alignment_length == 0:
         print("Warning: MSA has zero length. Cannot calculate diversity.")
         return {
            "column_hill": [],
            "mean_diversity": None,
            "normalized_diversity": None
        }

    for i in range(alignment_length):
        # Extract the i-th column as a string
        column_str = msa[:, i] # This slices the alignment to get the i-th column
        
        # Count character frequencies, ignoring gaps and stop codons, converting to uppercase
        chars_in_column = [char for char in column_str.upper() if char not in ['-', '.', '*']]

        if not chars_in_column:
            hill_div_col = 1.0 # Conserved 'non-information'
        else:
            freqs_counter = Counter(chars_in_column)
            total_chars = len(chars_in_column)
            freqs = [count / total_chars for count in freqs_counter.values()]

            entropy_col = -sum(p * math.log(p) for p in freqs if p > 0)
            hill_div_col = math.exp(entropy_col)

        column_hill.append(hill_div_col)

    if column_hill:
        mean_div = np.mean(column_hill)
    else:
        mean_div = None

    if mean_div is not None and S > 1:
         normalized_div = (mean_div - 1) / (S - 1)
         normalized_div = max(0.0, min(1.0, normalized_div))
    else:
         normalized_div = None

    return {
        "column_hill": column_hill,
        "mean_diversity": mean_div,
        "normalized_diversity": normalized_div
    }

# --- Example Usage: Process files in directories ---
if __name__ == "__main__":
    # Define input directories for trimmed alignments and trees
    # ADJUST THESE PATHS TO MATCH YOUR DIRECTORY STRUCTURE
    trimmed_alignments_dir = "intra_og_analysis/trimal_output_intra_og"
    fasttree_trees_dir = "intra_og_analysis/fasttree_output_intra_og_trees" # Assuming this is where the .nwk files are saved

    # Output file for results
    output_results_csv = "orthogroup_diversity_metrics.csv"

    # List to store results for all orthogroups
    all_og_diversity_results = []

    print(f"Processing trimmed alignments from: {trimmed_alignments_dir}")
    print(f"Processing phylogenetic trees from: {fasttree_trees_dir}")

    # Find all trimmed alignment files (assuming .fasta extension from trimal script)
    alignment_files = glob.glob(os.path.join(trimmed_alignments_dir, "*_trimmed.fasta"))
    print(f"Found {len(alignment_files)} trimmed alignment files.")

    if not alignment_files:
        print("No trimmed alignment files found. Please check the input directory.")
    else:
        # Process each alignment and its corresponding tree
        for msa_filepath in alignment_files:
            og_id = os.path.basename(msa_filepath).replace("_trimmed.fasta", "")
            tree_filepath = os.path.join(fasttree_trees_dir, f"{og_id}_fasttree.nwk") # Assuming .nwk extension from fasttree script

            print(f"\nProcessing Orthogroup: {og_id}")

            # --- Calculate MSA Diversity ---
            msa_diversity_results = {
                "OG_ID": og_id,
                "MSA_mean_diversity": None,
                "MSA_normalized_diversity": None,
                # "MSA_column_hill": None # Keeping this commented out to avoid large memory usage
            }
            try:
                msa = AlignIO.read(msa_filepath, "fasta")
                msa_results = hill_diversity_msa_py(msa)
                msa_diversity_results["MSA_mean_diversity"] = msa_results["mean_diversity"]
                msa_diversity_results["MSA_normalized_diversity"] = msa_results["normalized_diversity"]
                # Decide if you want to store the per-column diversity list.
                # It can be very large if you have many long alignments.
                # msa_diversity_results["MSA_column_hill"] = msa_results["column_hill"]
                if msa_results["mean_diversity"] is not None:
                     print(f"  MSA Diversity: Mean={msa_results['mean_diversity']:.4f}, Normalized={msa_results['normalized_diversity']:.4f}")
                else:
                     print("  MSA Diversity: Could not calculate.")
            except FileNotFoundError:
                print(f"  MSA file not found: {msa_filepath}")
            except Exception as e:
                print(f"  Error processing MSA file {msa_filepath}: {e}")
                # print(f"  Detailed error: {e}") # Uncomment for more detailed error info

            # --- Calculate Tree Diversity ---
            tree_diversity_results = {
                "Tree_Hill_diversity": None,
                "Tree_Norm_by_tips": None,
                "Tree_Norm_by_PD": None,
                "Tree_n_tips": None,
                "Tree_entropy": None,
            }
            if os.path.exists(tree_filepath): # Check if tree file exists (FastTree might fail for some OGs)
                try:
                    tree = Phylo.read(tree_filepath, "newick")
                    tree_results = hill_diversity_tree_py(tree)
                    tree_diversity_results["Tree_Hill_diversity"] = tree_results["Hill_diversity"]
                    tree_diversity_results["Tree_Norm_by_tips"] = tree_results["Norm_by_tips"]
                    tree_diversity_results["Tree_Norm_by_PD"] = tree_results["Norm_by_PD"]
                    tree_diversity_results["Tree_n_tips"] = tree_results["n_tips"]
                    tree_diversity_results["Tree_entropy"] = tree_results["entropy"]
                    if tree_results["Hill_diversity"] is not None:
                         print(f"  Tree Diversity: Hill={tree_results['Hill_diversity']:.4f}, Norm by Tips={tree_results['Norm_by_tips']:.4f}, Norm by PD={tree_results['Norm_by_PD']:.4f}, Tips={tree_results['n_tips']}")
                    else:
                         print("  Tree Diversity: Could not calculate.")
                except FileNotFoundError:
                     print(f"  Tree file not found (should not happen if check passed): {tree_filepath}")
                except Exception as e:
                    print(f"  Error processing tree file {tree_filepath}: {e}")
                    # print(f"  Detailed error: {e}") # Uncomment for more detailed error info
            else:
                 print(f"  Tree file not found for OG {og_id}: {tree_filepath}. Skipping tree diversity calculation.")


            # Combine results for this orthogroup
            og_results = {"OG_ID": og_id}
            og_results.update(msa_diversity_results)
            og_results.update(tree_diversity_results)

            all_og_diversity_results.append(og_results)

        # --- Store Results in a DataFrame and Save ---
        if all_og_diversity_results:
            df_diversity_results = pd.DataFrame(all_og_diversity_results)

            # Reorder columns for clarity (optional)
            ordered_cols = [
                "OG_ID", "Tree_n_tips", "Tree_Hill_diversity", "Tree_Norm_by_tips", "Tree_Norm_by_PD", "Tree_entropy",
                "MSA_mean_diversity", "MSA_normalized_diversity"
                # Add "MSA_column_hill" here if you chose to store it
            ]
            # Filter for columns that actually exist in the DataFrame
            ordered_cols_present = [col for col in ordered_cols if col in df_diversity_results.columns]
            df_diversity_results = df_diversity_results[ordered_cols_present]


            print(f"\nSaving diversity results for {len(df_diversity_results)} orthogroups to '{output_results_csv}'")
            try:
                df_diversity_results.to_csv(output_results_csv, index=False)
                print("Successfully saved diversity results.")
            except Exception as e:
                print(f"ERROR: Failed to save diversity results to {output_results_csv}: {e}")
        else:
            print("\nNo diversity results were calculated for any orthogroup.")

    print("\n--- Diversity Analysis Script Finished ---")
