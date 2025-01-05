import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import numpy as np
import os
import ast
import seaborn as sns


def load_results(filepath):
    """Load experiment results from TSV file"""
    return pd.read_csv(filepath, sep='\t')


def create_line_plots(df, output_dir):
    metrics = ['N50', 'Contig Count', 'Largest Contig']
    params = ['n', 'l', 'p']

    for metric in metrics:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle(f'{metric} vs Parameters')

        for idx, param in enumerate(params):
            # Create separate lines for other parameters combinations
            other_params = [p for p in params if p != param]

            for val1 in df[other_params[0]].unique():
                for val2 in df[other_params[1]].unique():
                    mask = (df[other_params[0]] == val1) & (df[other_params[1]] == val2)
                    data = df[mask].sort_values(param)
                    axes[idx].plot(data[param], data[metric],
                                   label=f'{other_params[0]}={val1}, {other_params[1]}={val2}')

            axes[idx].set_xlabel(param)
            axes[idx].set_ylabel(metric)
            axes[idx].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            axes[idx].grid(True)

        plt.tight_layout()
        plt.savefig(Path(output_dir) / f'{metric.lower().replace(" ", "_")}_analysis.png',
                    bbox_inches='tight', dpi=300)
        plt.close()


def visualize_results(results_file, contig_lengths_file, output_dir):
    """Main function to create all visualizations"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    df = load_results(results_file)

    create_line_plots(df, output_dir)
    # plot_contig_length_distribution(contig_lengths_file, output_dir)
    plot_contig_length_distribution2(contig_lengths_file, output_dir)

def plot_contig_length_distribution(contig_lengths_file, output_dir):
    df = pd.read_csv(contig_lengths_file, sep='\t')

    plt.figure(figsize=(12, 6))

    p_values = df['p'].unique()
    colors = plt.cm.viridis(np.linspace(0, 1, len(p_values)))

    max_length = max(max(ast.literal_eval(lengths)) for lengths in df['Contig Lengths'])
    bins = np.logspace(np.log10(10), np.log10(max_length), 50)

    for p, color in zip(p_values, colors):
        lengths = ast.literal_eval(df[df['p'] == p]['Contig Lengths'].iloc[0])
        plt.hist(lengths, bins=bins, alpha=0.7, label=f"p={p}",
                 density=True, color=color, edgecolor="black")

    plt.xscale('log')
    plt.xlabel("Contig Length (log scale)")
    plt.ylabel("Density")
    plt.title("Contig Length Distribution for Different Error Probabilities")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    output_path = os.path.join(output_dir, "contig_length_distribution.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Contig length distribution plot saved to '{output_path}'")


def plot_contig_length_distribution2(contig_lengths_file, output_dir):
    # Read the data
    df = pd.read_csv(contig_lengths_file, sep='\t')

    # Convert the string representation of lists to actual lists
    df['Contig Lengths'] = df['Contig Lengths'].apply(ast.literal_eval)

    # Create a new DataFrame with exploded lengths
    df_exploded = df.explode('Contig Lengths').reset_index(drop=True)
    df_exploded['Contig Lengths'] = df_exploded['Contig Lengths'].astype(float)

    # Create the plot
    plt.figure(figsize=(12, 6))
    sns.kdeplot(data=df_exploded, x='Contig Lengths', hue='p', fill=True, alpha=0.4)

    plt.title('Contig Length Distribution by P-value', fontsize=14)
    plt.xlabel('Contig Length', fontsize=12)
    plt.ylabel('Density', fontsize=12)

    # Save the plot
    plt.savefig(os.path.join(output_dir, 'contig_length_distribution.png'))
    plt.close()