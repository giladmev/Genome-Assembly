from collections import defaultdict
from Bio import SeqIO


def build_overlap_graph(reads, min_overlap):
    graph = defaultdict(list)
    for i, read1 in enumerate(reads):
        for j, read2 in enumerate(reads):
            if i != j:
                overlap = find_overlap(read1, read2, min_overlap)
                if overlap >= min_overlap:
                    graph[read1].append((read2, overlap))
    return graph


def find_overlap(read1, read2, min_overlap):
    for i in range(len(read1) - min_overlap + 1):
        if read2.startswith(read1[i:]):
            return len(read1) - i
    return 0


def assemble_genome(graph):
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
