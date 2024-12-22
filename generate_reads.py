import random
from Bio import SeqIO

def load_genome_sequence(fasta_file):
    record = next(SeqIO.parse(fasta_file, "fasta"))
    genome_sequence = str(record.seq)
    genome_length = len(genome_sequence)
    return genome_sequence, genome_length

def save_reads_to_fasta(reads, output_file_name):
    with open(output_file_name, "w") as output_file:
        for i, read in enumerate(reads):
            output_file.write(f">read_{i + 1}\n{read}\n")
    print(f"Generated reads saved to '{output_file_name}'\n")

def generate_reads(fasta_file, read_length, num_reads, error_prob=0):
    print(f"Generating reads from '{fasta_file}'")
    genome_sequence, genome_length = load_genome_sequence(fasta_file)
    coverage = (num_reads * read_length) / genome_length
    print(f"Genome Length: {genome_length} bp")
    print(f"Desired Number of Reads: {num_reads}")
    print(f"Read Length: {read_length} bp")
    print(f"Error Probability per Base: {error_prob}")
    print(f"Approximate Coverage: {coverage}x")

    def introduce_error(base):
        bases = ["A", "T", "C", "G"]
        bases.remove(base)
        return random.choice(bases)

    reads = []
    for _ in range(num_reads):
        start_pos = random.randint(0, genome_length - read_length)
        read = genome_sequence[start_pos:start_pos + read_length]

        if error_prob > 0:
            read_with_errors = []
            for base in read:
                if random.random() < error_prob:
                    read_with_errors.append(introduce_error(base))
                else:
                    read_with_errors.append(base)
            read = "".join(read_with_errors)

        reads.append(read)
    base_path = "data\\output\\"
    output_file_name = f"{base_path}reads_l{read_length}_n{num_reads}_p{int(error_prob*100)}%.fasta"
    save_reads_to_fasta(reads, output_file_name)
    return output_file_name

