import itertools
from app.generate_reads import generate_reads
from app.assamble_genome_improving_performance import assemble_reads
from app.calculate_assembled_metrics import calculate_n50, calculate_largest_contig_fraction
from app.visualize_results import visualize_results
import yaml
from pathlib import Path
import os
from Bio import SeqIO


def run_all_experiments(num_reads_values, read_length_values, error_prob_values, genome_fasta_file, output_dir):
    results = []
    for n, l, p in itertools.product(num_reads_values, read_length_values, error_prob_values):
        is_l_n_max = n == max(num_reads_values) and l == max(read_length_values)
        result = run_experiment(n, l, p, genome_fasta_file, output_dir, document_contig_lengths=is_l_n_max)
        results.append(result)
    return results


def run_experiment(n, l, p, genome_fasta_file, output_dir, document_contig_lengths=False):
    # Generate reads with specified parameters
    reads_file_path = generate_reads(genome_fasta_file, l, n, output_dir, error_prob=p)

    # Assemble reads
    assembled_contigs_file_path = assemble_reads(reads_file_path, output_dir, min_overlap=20)
    contigs = list(SeqIO.parse(assembled_contigs_file_path, "fasta"))

    n50 = calculate_n50(contigs)
    contig_count = len(contigs)
    largest_contig = max(len(contig) for contig in contigs)
    genome_fraction = calculate_largest_contig_fraction(contigs, genome_fasta_file)

    results = {
        "n": n,
        "l": l,
        "p": p,
        "N50": n50,
        "Contig Count": contig_count,
        "Largest Contig": largest_contig,
        "Genome Fraction": genome_fraction
    }

    if document_contig_lengths:
        contig_lengths = [len(contig) for contig in contigs]
        results["Contig Lengths"] = contig_lengths

    return results


def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    # Create output directory if it doesn't exist
    Path(config['output']['directories']['faste_files']).mkdir(parents=True, exist_ok=True)
    Path(config['output']['directories']['results']).mkdir(parents=True, exist_ok=True)
    return config


def document_findings(results, results_output_dir):
    results_file = os.path.join(results_output_dir, "experiment_results.txt")
    with open(results_file, "w") as f:
        f.write("n\tl\tp\tN50\tContig Count\tLargest Contig\tGenome Fraction\n")
        for result in results:
            f.write(
                f"{result['n']}\t{result['l']}\t{result['p']}\t{result['N50']}\t{result['Contig Count']}\t{result['Largest Contig']}\t{result['Genome Fraction']}\n")

    print("Experiment results saved to 'experiment_results.txt'")
    return results_file

def document_contig_lengths(results, results_output_dir):
    contig_lengths_file = os.path.join(results_output_dir, "contig_lengths.txt")
    with open(contig_lengths_file, "w") as f:
        f.write("p\tContig Lengths\n")
        for result in results:
            p = result['p']
            if 'Contig Lengths' in result:
                contig_lengths = result['Contig Lengths']
                f.write(f"{p}\t{contig_lengths}\n")
    print("Contig lengths saved to 'contig_lengths.txt'")
    return contig_lengths_file



if __name__ == "__main__":
    config = load_config()

    genome_fasta_file_path = config['input']['genome_fasta']
    fasta_files_output_dir = config['output']['directories']['faste_files']
    results_dir = config['output']['directories']['results']
    n_values = config['sequencing']['num_reads']
    l_values = config['sequencing']['read_length']
    p_values = config['assembly']['error_probabilities']

    # Run experiments for all combinations of parameters
    experiments_results = run_all_experiments(n_values, l_values, p_values, genome_fasta_file_path, fasta_files_output_dir)

    results_file_path = document_findings(experiments_results, results_dir)
    contig_lengths_file_path = document_contig_lengths(experiments_results, results_dir)

    contig_lengths_file_path = os.path.join(results_dir, "contig_lengths.txt")
    results_file_path = os.path.join(results_dir, "experiment_results.txt")
    visualize_results(results_file_path, contig_lengths_file_path, results_dir)
