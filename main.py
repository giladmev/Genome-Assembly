import time
from generate_reads import generate_reads
from assamble_genome_improving_performance import assemble_reads
from generate_assembly_report import generate_assembly_report
from filtering import filter_contigs_parallel


if __name__ == "__main__":
    genome_fasta_file = "data/input/phiX_genome.fasta"
    read_length = 100
    num_reads = 10000
    min_overlap = 20
    error_prob = 0.01

    ASSAMBLE_CONTIGS = True
    if ASSAMBLE_CONTIGS:
        # Generate error-free reads
        error_free_reads_file = generate_reads(genome_fasta_file, read_length, num_reads)

        # Generate error-prone reads
        error_prone_reads_file = generate_reads(genome_fasta_file, read_length, num_reads, error_prob=error_prob)

        # Assemble the reads
        start_time = time.time()
        contigs_error_free = assemble_reads(error_free_reads_file, min_overlap)
        filter_contigs_parallel(contigs_error_free)
        print(f"Time taken to assemble error-free reads: {time.time() - start_time:.2f} seconds\n")

        start_time = time.time()
        contigs_error_prone = assemble_reads(error_prone_reads_file, min_overlap)
        filter_contigs_parallel(contigs_error_prone)
        print(f"Time taken to assemble error-prone reads: {time.time() - start_time:.2f} seconds\n")
    else:
        contigs_error_free = "data/output/assembled_contigs_l100_n10000_p0%.fasta"
        contigs_error_prone = "data/output/assembled_contigs_l100_n10000_p1%.fasta"

    # Generate assembly report
    generate_assembly_report(contigs_error_free, contigs_error_prone, genome_fasta_file)
