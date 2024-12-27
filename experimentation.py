import itertools
from app.generate_reads import generate_reads
from app.assamble_genome_improving_performance import assemble_reads
from app.generate_assembly_report import calculate_n50
import yaml
from pathlib import Path
import os
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO


def run_experiment(n, l, p, genome_fasta_file, output_dir):
    # Generate reads with specified parameters
    reads_file_path = generate_reads(genome_fasta_file, l, n, output_dir, error_prob=p)

    file_details = '_'.join(reads_file_path.split(os.path.sep)[-1].split(".")[0].split("_")[1:])
    file_name = f"assembled_contigs_{file_details}.fasta"
    assembled_contigs_file_path = os.path.join(output_dir, file_name)

    if not os.path.exists(assembled_contigs_file_path):
        # Assemble reads and compute metrics
        assembled_contigs_file_path = assemble_reads(reads_file_path, output_dir, min_overlap=20)

    # filter_contigs(assembled_contigs_file_path)
    contigs = list(SeqIO.parse(assembled_contigs_file_path, "fasta"))

    n50 = calculate_n50(contigs)
    contig_count = len(contigs)
    largest_contig = max(len(contig) for contig in contigs)

    # Print results
    print(f"Parameters: n={n}, l={l}, p={p}")
    print(f"N50: {n50}, Contig Count: {contig_count}, Largest Contig: {largest_contig}\n")

    return {
        "n": n,
        "l": l,
        "p": p,
        "N50": n50,
        "Contig Count": contig_count,
        "Largest Contig": largest_contig
    }

def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    # Create output directory if it doesn't exist
    Path(config['output']['directories']['faste_files']).mkdir(parents=True, exist_ok=True)
    Path(config['output']['directories']['results']).mkdir(parents=True, exist_ok=True)
    return config

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

if __name__ == "__main__":
    config = load_config()

    genome_fasta_file = config['input']['genome_fasta']
    fasta_files_output_dir = config['output']['directories']['faste_files']
    results_output_dir = config['output']['directories']['results']
    n_values = config['sequencing']['num_reads']
    l_values = config['sequencing']['read_length']
    p_values = config['assembly']['error_prob']

    # Run experiments for all combinations of parameters
    results_file = os.path.join(results_output_dir, "experiment_results.txt")
    results = []

    if os.path.exists(results_file):
        print(f"Results file '{results_file}' found. Reading results...")
        visualize_results(results_file, results_output_dir)
    else:
        for n, l, p in itertools.product(n_values, l_values, p_values):
            result = run_experiment(n, l, p, genome_fasta_file, fasta_files_output_dir)
            results.append(result)

        # Document findings
        with open(results_file, "w") as f:
            f.write("n\tl\tp\tN50\tContig Count\tLargest Contig\n")
            for result in results:
                f.write(
                    f"{result['n']}\t{result['l']}\t{result['p']}\t{result['N50']}\t{result['Contig Count']}\t{result['Largest Contig']}\n")

        print("Experiment results saved to 'experiment_results.txt'")
        visualize_results(results_file, results_output_dir)
