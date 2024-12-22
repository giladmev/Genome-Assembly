import time
from app.generate_reads import generate_reads
from app.assamble_genome_improving_performance import assemble_reads
from app.generate_assembly_report import generate_assembly_report
from app.filtering import filter_contigs
import yaml
from pathlib import Path


def load_config(config_path="config.yaml"):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    # Create output directory if it doesn't exist
    Path(config['output']['directory']).mkdir(parents=True, exist_ok=True)

    return config

if __name__ == "__main__":
    config = load_config()

    genome_fasta_file = config['input']['genome_fasta']
    output_dir = config['output']['directory']
    read_length = config['sequencing']['read_length']
    num_reads = config['sequencing']['num_reads']
    min_overlap = config['assembly']['min_overlap']
    error_prob = config['assembly']['error_prob']

    # Generate error-free reads
    error_free_reads_file = generate_reads(genome_fasta_file, read_length, output_dir, num_reads)

    # Generate error-prone reads
    error_prone_reads_file = generate_reads(genome_fasta_file, read_length, num_reads, output_dir, error_prob=error_prob)

    # Assemble the reads
    start_time = time.time()
    contigs_error_free = assemble_reads(error_free_reads_file, min_overlap, output_dir)
    filter_contigs(contigs_error_free)
    print(f"Time taken to assemble error-free reads: {time.time() - start_time:.2f} seconds\n")

    start_time = time.time()
    contigs_error_prone = assemble_reads(error_prone_reads_file, min_overlap, output_dir)
    filter_contigs(contigs_error_prone)
    print(f"Time taken to assemble error-prone reads: {time.time() - start_time:.2f} seconds\n")

    # Generate assembly report
    generate_assembly_report(contigs_error_free, contigs_error_prone, genome_fasta_file)
