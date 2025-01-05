
from Bio import SeqIO

def calculate_n50(contigs):
    total_length = sum(len(contig) for contig in contigs)
    sorted_contigs = sorted(contigs, key=len, reverse=True)
    cumulative_length = 0
    for contig in sorted_contigs:
        cumulative_length += len(contig)
        if cumulative_length >= total_length / 2:
            return len(contig)
    return 0


def calculate_genome_fraction(contigs, reference):
    assembly_seq = ''.join([str(contig.seq) for contig in contigs])
    covered_bases = 0
    for i in range(len(reference)):
        if i < len(assembly_seq) and reference[i] == assembly_seq[i]:
            covered_bases += 1
    return (covered_bases / len(reference)) * 100.


def calculate_largest_contig_fraction(contigs, reference_file):
    # Find the largest contig
    largest_contig = max((str(contig.seq) for contig in contigs), key=len)

    # Read the reference genome
    reference = str(next(SeqIO.parse(reference_file, "fasta")).seq)

    best_match_start = -1
    best_match_score = 0

    # Find the best matching position for the largest contig
    for i in range(len(reference) - len(largest_contig) + 1):
        match_score = sum(1 for a, b in zip(largest_contig, reference[i:]) if a == b)
        if match_score > best_match_score:
            best_match_score = match_score
            best_match_start = i

    # Calculate the number of matching bases
    if best_match_start != -1:
        matching_bases = sum(1 for a, b in zip(largest_contig, reference[best_match_start:]) if a == b)
    else:
        matching_bases = 0

    # Calculate the fraction
    largest_contig_fraction = (matching_bases / len(reference)) * 100

    return largest_contig_fraction
