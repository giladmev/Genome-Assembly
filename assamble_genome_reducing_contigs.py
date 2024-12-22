import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from itertools import combinations


def find_overlap(read1, read2, min_overlap, max_error_rate=0.05):
    for i in range(min(len(read1), len(read2)), min_overlap - 1, -1):
        mismatches = sum(a != b for a, b in zip(read1[-i:], read2[:i]))
        if mismatches / i <= max_error_rate:
            return i
    return 0


def find_overlaps_chunk(chunk, all_reads, min_overlap, max_error_rate):
    overlaps = []
    for read1 in chunk:
        for read2 in all_reads:
            if read1 != read2:
                overlap = find_overlap(read1, read2, min_overlap, max_error_rate)
                if overlap >= min_overlap:
                    overlaps.append((read1, read2, overlap))
    return overlaps


def find_overlaps_parallel(reads, min_overlap, max_error_rate=0.05):
    num_processes = mp.cpu_count()
    chunk_size = max(1, len(reads) // num_processes)
    chunks = [reads[i:i + chunk_size] for i in range(0, len(reads), chunk_size)]

    with mp.Pool(processes=num_processes) as pool:
        all_overlaps = pool.starmap(find_overlaps_chunk,
                                    [(chunk, reads, min_overlap, max_error_rate) for chunk in chunks])

    return [overlap for sublist in all_overlaps for overlap in sublist]


def build_overlap_graph(overlaps):
    graph = defaultdict(list)
    for read1, read2, overlap in overlaps:
        graph[read1].append((read2, overlap))
    return graph


def extend_contig(graph, start_read):
    contig = start_read
    current_read = start_read
    path = [start_read]
    while True:
        next_reads = [(read, overlap) for read, overlap in graph[current_read] if read not in path]
        if not next_reads:
            break
        next_read, overlap = max(next_reads, key=lambda x: x[1])
        contig += next_read[overlap:]
        path.append(next_read)
        current_read = next_read
    return contig


def assemble_genome_parallel(graph, min_contig_length):
    with mp.Pool(processes=mp.cpu_count()) as pool:
        contigs = pool.starmap(extend_contig, [(graph, start_read, min_contig_length) for start_read in graph.keys()])
    return [contig for contig in contigs if contig]


def filter_contigs(contigs, min_length=100, similarity_threshold=0.9):
    filtered_contigs = [contig for contig in contigs if len(contig) >= min_length]
    filtered_contigs.sort(key=len, reverse=True)

    final_contigs = []
    for contig in filtered_contigs:
        if not any(is_similar(contig, fc, similarity_threshold) for fc in final_contigs):
            final_contigs.append(contig)

    return final_contigs


def is_similar(contig1, contig2, threshold):
    shorter, longer = sorted([contig1, contig2], key=len)
    for i in range(len(longer) - len(shorter) + 1):
        if sum(a == b for a, b in zip(shorter, longer[i:i + len(shorter)])) / len(shorter) >= threshold:
            return True
    return False


def assemble_reads(fasta_file, min_overlap, min_contig_length=100, max_error_rate=0.05):
    reads = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    print("Assembling reads...")
    print(f"Number of reads: {len(reads)}")
    print(f"Minimum overlap: {min_overlap}")

    overlaps = find_overlaps_parallel(reads, min_overlap, max_error_rate)
    graph = build_overlap_graph(overlaps)
    contigs = assemble_genome_parallel(graph, min_contig_length)
    filtered_contigs = filter_contigs(contigs, min_contig_length)

    print(f"Number of contigs: {len(filtered_contigs)}")
    print(f"Longest contig length: {max(len(contig) for contig in filtered_contigs)}")

    file_details = '_'.join(fasta_file.split("\\")[-1].split(".")[0].split("_")[1:])
    assembled_contigs_file = f"data\\output\\assembled_contigs_{file_details}.fasta"
    with open(assembled_contigs_file, "w") as f:
        for i, contig in enumerate(filtered_contigs):
            f.write(f">contig_{i + 1}\n{contig}\n")

    print(f"Assembled contigs saved to '{assembled_contigs_file}'")
    return assembled_contigs_file
