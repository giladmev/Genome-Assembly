import time
from app.generate_reads import generate_reads
from app.assamble_genome_improving_performance import assemble_reads
from app.calculate_assembled_metrics import generate_assembly_report
import yaml
from pathlib import Path


def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    Path(config['output']['directories']['faste_files']).mkdir(parents=True, exist_ok=True)
    Path(config['output']['directories']['results']).mkdir(parents=True, exist_ok=True)

    return config

if __name__ == "__main__":
    config = load_config()

    genome_fasta_file = config['input']['genome_fasta']
    fasta_files_output_dir = config['output']['directories']['faste_files']
    results_output_dir = config['output']['directories']['results']
    read_length = config['sequencing']['read_length'][1]
    num_reads = config['sequencing']['num_reads'][1]
    min_overlap = config['assembly']['min_overlap']
    error_probabilities = config['assembly']['error_probabilities'][1]

    # Generate error-free reads
    error_free_reads_file_path = generate_reads(genome_fasta_file, read_length, num_reads, fasta_files_output_dir)

    # Generate error-prone reads
    error_prone_reads_file_path = generate_reads(genome_fasta_file, read_length, num_reads, fasta_files_output_dir, error_prob=error_prob)

    # Assemble the reads
    start_time = time.time()
    contigs_error_free_file_path = assemble_reads(error_free_reads_file_path, fasta_files_output_dir, min_overlap)
    print(f"Time taken to assemble error-free reads: {time.time() - start_time:.2f} seconds\n")

    start_time = time.time()
    contigs_error_prone_file_path = assemble_reads(error_prone_reads_file_path, fasta_files_output_dir, min_overlap)
    print(f"Time taken to assemble error-prone reads: {time.time() - start_time:.2f} seconds\n")

    # Generate assembly report
    generate_assembly_report(contigs_error_free_file_path, contigs_error_prone_file_path, genome_fasta_file, results_output_dir)
