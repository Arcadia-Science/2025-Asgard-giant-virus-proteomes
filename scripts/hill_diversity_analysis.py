import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse  # ADDED: Import argparse

def hill_diversity(p, q):
    """Calculates Hill diversity."""
    p = p[p > 0]
    if q == 1:
        return np.exp(-np.sum(p * np.log(p)))
    else:
        return np.sum(p**q)**(1/(1-q))

def main():
    # ADDED: Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate and plot Hill diversity from orthogroup data.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input orthogroup diversity metrics CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output plot PNG file.")
    args = parser.parse_args()

    # MODIFIED: Use the input file from command-line arguments
    df = pd.read_csv(args.input)

    q_values = np.linspace(0, 3, 50)
    results = []

    for index, row in df.iterrows():
        proportions = row.values[1:] / np.sum(row.values[1:])
        for q in q_values:
            diversity = hill_diversity(proportions, q)
            results.append({
                'Orthogroup': row['Orthogroup'],
                'q': q,
                'Diversity': diversity
            })

    results_df = pd.DataFrame(results)

    # Plotting
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))

    sns.lineplot(data=results_df, x='q', y='Diversity', hue='Orthogroup', ax=ax, palette="viridis")

    ax.set_title('Hill Diversity of Eukaryotic Domains in Orthogroups', fontsize=16)
    ax.set_xlabel('Order q', fontsize=12)
    ax.set_ylabel('Diversity (Effective Number of Domains)', fontsize=12)
    ax.legend(title='Orthogroup', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    # MODIFIED: Save the plot to the specified output file
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to {args.output}")

if __name__ == '__main__':
    main()