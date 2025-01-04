import itertools
from app.generate_reads import generate_reads
from app.assamble_genome_improving_performance import assemble_reads
from app.generate_assembly_report import calculate_n50
from app.visualize_results import visualize_results
import yaml
from pathlib import Path
import os
from Bio import SeqIO


def run_all_experiments(num_reads_values, read_length_values, error_prob_values, genome_fasta_file, output_dir):
    results = []
    for n, l, p in itertools.product(num_reads_values, read_length_values, error_prob_values):
        result = run_experiment(n, l, p, genome_fasta_file, output_dir)
        results.append(result)
    return results

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

def document_findings(results):
    results_file = os.path.join(results_output_dir, "experiment_results.txt")
    with open(results_file, "w") as f:
        f.write("n\tl\tp\tN50\tContig Count\tLargest Contig\n")
    for result in results:
        f.write(f"{result['n']}\t{result['l']}\t{result['p']}\t{result['N50']}\t{result['Contig Count']}\t{result['Largest Contig']}\n")

    print("Experiment results saved to 'experiment_results.txt'")
    return results_file


if __name__ == "__main__":
    config = load_config()

    genome_fasta_file = config['input']['genome_fasta']
    fasta_files_output_dir = config['output']['directories']['faste_files']
    results_output_dir = config['output']['directories']['results']
    n_values = config['sequencing']['num_reads']
    l_values = config['sequencing']['read_length']
    p_values = config['assembly']['error_prob']

    # Run experiments for all combinations of parameters
    experiments_results = run_all_experiments(n_values, l_values, p_values, genome_fasta_file, fasta_files_output_dir)

    results_file_path = document_findings(experiments_results)
    visualize_results(results_file_path, results_output_dir)
