from generate_reads import generate_reads
from assamble_genome import assemble_reads
from generate_assembly_report import generate_assembly_report


if __name__ == "__main__":
    # fasta_file = "phiX_genome.fasta"
    # read_length = 100
    # num_reads = 10000
    # min_overlap = 20
    # error_prob = 0.01
    #
    # # Generate error-free reads
    # error_free_reads_file = generate_reads(fasta_file, read_length, num_reads)
    #
    # # Generate error-prone reads
    # error_prone_reads_file = generate_reads(fasta_file, read_length, num_reads, error_prob=error_prob)
    #
    # # Assemble the reads
    # contigs_error_free = assemble_reads(error_free_reads_file, min_overlap)
    # contigs_error_prone = assemble_reads(error_prone_reads_file, min_overlap)
    generate_assembly_report("assembled_contigs_error_free_reads.fasta", "assembled_contigs_error_prone_reads.fasta", "phiX_genome.fasta")