import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


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

def visualize_results(results_file, output_dir):
    """Main function to create all visualizations"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    df = load_results(results_file)

    create_line_plots(df, output_dir)