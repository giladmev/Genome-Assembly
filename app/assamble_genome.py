from collections import defaultdict
from Bio import SeqIO


def build_overlap_graph(reads, min_overlap):
    """
     Constructs an overlap graph where nodes are reads and edges represent overlaps between reads.

    Args:
        reads (list): List of read sequences (str).
        min_overlap (int): Minimum overlap length to consider.

    Returns:
        dict: A dictionary representing the overlap graph.
    """
    graph = defaultdict(list)
    for i, read1 in enumerate(reads):
        for j, read2 in enumerate(reads):
            if i != j:
                overlap = find_overlap(read1, read2, min_overlap)
                if overlap >= min_overlap:
                    graph[read1].append((read2, overlap))
    return graph


def find_overlap(read1, read2, min_overlap):
    """
    Finds the overlap between two reads..

    Args:
        read1 (str): The first read sequence.
        read2 (str): The second read sequence.
        min_overlap (int): Minimum overlap length to consider.

    Returns:
        int: The length of the overlap.
    """
    for i in range(len(read1) - min_overlap + 1):
        if read2.startswith(read1[i:]):
            return len(read1) - i
    return 0


def assemble_genome(graph):
    """
    Traverses the graph to generate contigs.

    Args:
        graph (dict): The overlap graph.

    Returns:
        list: A list of assembled contigs (str).
    """
    visited = set()
    contigs = []

    for start_read in list(graph.keys()):
        if start_read not in visited:
            contig = [start_read]
            visited.add(start_read)
            current_read = start_read

            while True:
                next_reads = graph[current_read]
                if not next_reads:
                    break

                next_read, overlap = max(next_reads, key=lambda x: x[1])
                if next_read in visited:
                    break

                contig.append(next_read[overlap:])
                visited.add(next_read)
                current_read = next_read

            contigs.append("".join(contig))

    return contigs

def assemble_reads(fasta_file, min_overlap):
    """
    Main function that loads reads, builds the graph, and assembles the genome.

    Args:
        fasta_file (str): Path to the FASTA file containing the reads.
        min_overlap (int): Minimum overlap length to consider.

    Returns:
        str: Path to the output FASTA file containing the assembled contigs.
    """
    print(f"Assembling reads from '{fasta_file}'")
    reads = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

    print(f"Number of reads: {len(reads)}")
    print(f"Minimum overlap: {min_overlap}")

    graph = build_overlap_graph(reads, min_overlap)
    contigs = assemble_genome(graph)

    print(f"Number of contigs: {len(contigs)}")
    print(f"Longest contig length: {max(len(contig) for contig in contigs)}")

    # Save contigs to a file
    assembled_contigs_file = f"assembled_contigs_{fasta_file}"
    with open(assembled_contigs_file, "w") as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i + 1}\n{contig}\n")
    print(f"Assembled contigs saved to '{assembled_contigs_file}'")
    print()
    return assembled_contigs_file
