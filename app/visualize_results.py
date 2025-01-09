import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import numpy as np
import os
import ast
import seaborn as sns
from scipy import stats
from scipy.stats import gaussian_kde


def load_results(filepath):
    """Load experiment results from TSV file"""
    return pd.read_csv(filepath, sep='\t')

def visualize_results(results_file, contig_lengths_file, output_dir):
    """Main function to create all visualizations"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    plot_genome_fraction_vs_n_l(results_file, output_dir)
    plot_contig_length_distribution(contig_lengths_file, output_dir)
    plot_metrics_by_p(results_file, output_dir)

def plot_genome_fraction_vs_n_l(results_file, output_dir):
    df = load_results(results_file)
    # Calculate n * l
    df['n*l'] = df['n'] * df['l']

    # Plot the data
    plt.figure(figsize=(12, 6))

    # Create a line plot with different colors for each p value
    sns.lineplot(
        data=df,
        x='n*l',
        y='Genome Fraction',
        hue='p',
        palette='tab10',  # Choose a palette for colors
        marker='o',
        errorbar=None
    )

    # Customize the plot
    plt.title('Genome Fraction vs n*l (Colored by p)', fontsize=16)
    plt.xlabel('n * l', fontsize=12)
    plt.ylabel('Genome Fraction', fontsize=12)
    plt.legend(title='p', fontsize=10, title_fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)

    output_path = os.path.join(output_dir, "genome_fraction_by_nl.png")
    plt.savefig(output_path)
    plt.close()
    print(f"genome_fraction_by_n*l plot saved to '{output_path}'")

def plot_contig_length_distribution(contig_lengths_file, output_dir):
    # Read the data
    df = pd.read_csv(contig_lengths_file, sep='\t')

    # Convert the string representation of lists to actual lists
    df['Contig Lengths'] = df['Contig Lengths'].apply(ast.literal_eval)

    # Create a new DataFrame with exploded lengths
    df_exploded = df.explode('Contig Lengths').reset_index(drop=True)
    df_exploded['Contig Lengths'] = df_exploded['Contig Lengths'].astype(float)

    # Extract unique p-values
    unique_p_values = df['p'].unique()

    # Define the color palette
    palette = sns.color_palette("Blues", n_colors=len(unique_p_values))[::-1]
    palette[unique_p_values.tolist().index(unique_p_values[0])] = (0.0, 1.0, 0.0)  # Green for the first p value

    # Create the plot
    plt.figure(figsize=(12, 6))
    sns.kdeplot(data=df_exploded, x='Contig Lengths', hue='p', fill=True, alpha=0.4, palette=palette)
    plt.xscale('log')

    plt.xlim(1, 5500)
    plt.title('Contig Length Distribution by P-value', fontsize=14)
    plt.xlabel('Contig Length', fontsize=12)
    plt.ylabel('Density', fontsize=12)

    # Add labels at the peak density of each group
    for p_value, color in zip(unique_p_values, palette):
        # Filter the data for the current p value
        data = df_exploded[df_exploded['p'] == p_value]['Contig Lengths']

        # Compute the KDE using scipy
        kde = gaussian_kde(data)
        x_vals = np.linspace(0, 1200, 1000)  # Define range for x values
        y_vals = kde(x_vals)

        # Find the peak (highest density)
        max_idx = np.argmax(y_vals)
        x_peak, y_peak = x_vals[max_idx], y_vals[max_idx]

        # Annotate the plot
        plt.text(x_peak, y_peak, f'p={int(p_value)}', fontsize=10, ha='center', color=color)

    x_label_y_offset = plt.ylim()[0] - 0.02 * (plt.ylim()[1] - plt.ylim()[0])  # Offset below the x-axis
    plt.text(5500, x_label_y_offset, '5500', fontsize=10, ha='center', va='top')

    # Save the plot
    output_path = os.path.join(output_dir, 'contig_length_distribution.png')
    plt.savefig(output_path)
    print(f"contig_length_distribution plot saved to '{output_path}'")


    plt.close()

def plot_metrics_by_p(results_file, output_dir):
    df = load_results(results_file)
    # Find maximum n and l values
    max_n = df['n'].max()
    max_l = df['l'].max()

    # Filter data for max n and l
    df_filtered = df[(df['n'] == max_n) & (df['l'] == max_l)]

    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 6))

    # Colors for different p values
    colors = sns.color_palette("husl", len(df_filtered['p'].unique()))

    # Width of bars
    bar_width = 0.4

    # Create positions for bars
    p_values = sorted(df_filtered['p'].unique())
    x = np.arange(len(p_values))

    # Plot N50
    ax1.bar(x, df_filtered['N50'], bar_width, label='N50', color=colors)
    ax1.set_title('N50')
    ax1.set_xlabel('p-value')
    ax1.set_xticks(x)
    ax1.set_xticklabels(p_values)

    # Plot Largest Contig
    ax2.bar(x, df_filtered['Largest Contig'], bar_width, label='Largest Contig', color=colors)
    ax2.set_title('Largest Contig')
    ax2.set_xlabel('p-value')
    ax2.set_xticks(x)
    ax2.set_xticklabels(p_values)

    # Plot Genome Fraction
    ax3.bar(x, df_filtered['Genome Fraction'], bar_width, label='Genome Fraction', color=colors)
    ax3.set_title('Genome Fraction')
    ax3.set_xlabel('p-value')
    ax3.set_xticks(x)
    ax3.set_xticklabels(p_values)

    # Add main title
    plt.suptitle(f'Metrics Comparison for n={max_n}, l={max_l}', fontsize=14)

    # Adjust layout
    plt.tight_layout()

    output_path = os.path.join(output_dir, "metric_comparison_by_p.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"metric_comparison_by_p plot saved to '{output_path}'")
