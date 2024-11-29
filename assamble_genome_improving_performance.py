import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from itertools import combinations


def find_overlap(read1, read2, min_overlap):
    for i in range(min(len(read1), len(read2)), min_overlap - 1, -1):
        if read1.endswith(read2[:i]):
            return i
    return 0


def find_overlaps_chunk(chunk, all_reads, min_overlap):
    overlaps = []
    for read1 in chunk:
        for read2 in all_reads:
            if read1 != read2:
                overlap = find_overlap(read1, read2, min_overlap)
                if overlap >= min_overlap:
                    overlaps.append((read1, read2, overlap))
    return overlaps


def find_overlaps_parallel(reads, min_overlap):
    num_processes = mp.cpu_count()
    chunk_size = max(1, len(reads) // num_processes)
    chunks = [reads[i:i + chunk_size] for i in range(0, len(reads), chunk_size)]

    with mp.Pool(processes=num_processes) as pool:
        all_overlaps = pool.starmap(find_overlaps_chunk, [(chunk, reads, min_overlap) for chunk in chunks])

    return [overlap for sublist in all_overlaps for overlap in sublist]


def build_overlap_graph(overlaps):
    graph = defaultdict(list)
    for read1, read2, overlap in overlaps:
        graph[read1].append((read2, overlap))
    return graph


def assemble_genome_parallel(graph):
    def extend_contig(start_read):
        contig = [start_read]
        current_read = start_read
        while True:
            next_reads = graph[current_read]
            if not next_reads:
                break
            next_read, overlap = max(next_reads, key=lambda x: x[1])
            if next_read in contig:
                break
            contig.append(next_read[overlap:])
            current_read = next_read
        return "".join(contig)

    with mp.Pool(processes=mp.cpu_count()) as pool:
        contigs = pool.map(extend_contig, graph.keys())

    return [contig for contig in contigs if contig]


def assemble_reads(fasta_file, min_overlap):
    reads = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

    print(f"Number of reads: {len(reads)}")
    print(f"Minimum overlap: {min_overlap}")

    # Find overlaps in parallel
    overlaps = find_overlaps_parallel(reads, min_overlap)

    # Build overlap graph
    graph = build_overlap_graph(overlaps)

    # Assemble genome in parallel
    contigs = assemble_genome_parallel(graph)

    print(f"Number of contigs: {len(contigs)}")
    print(f"Longest contig length: {max(len(contig) for contig in contigs)}")

    assembled_contigs_file = "assembled_contigs_parallel.fasta"
    with open(assembled_contigs_file, "w") as f:
        for i, contig in enumerate(contigs):
            f.write(f">contig_{i + 1}\n{contig}\n")

    print(f"Assembled contigs saved to '{assembled_contigs_file}'")
    return contigs, assembled_contigs_file
