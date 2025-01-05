import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
import os
from app.filtering import filter_contigs


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

def extend_contig(graph, start_read):
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

def assemble_genome_parallel(graph):
    with mp.Pool(processes=mp.cpu_count()) as pool:
        contigs = pool.starmap(extend_contig, [(graph, start_read) for start_read in graph.keys()])

    return [contig for contig in contigs if contig]


def assemble_reads(fasta_file, output_dir, min_overlap, overwrite=False):
    file_details = '_'.join(fasta_file.split(os.path.sep)[-1].split(".")[0].split("_")[1:])
    file_name = f"assembled_contigs_{file_details}.fasta"
    assembled_contigs_file_path = os.path.join(output_dir, file_name)
    if not overwrite and os.path.exists(output_dir):
        output_files = os.listdir(output_dir)
        if file_name in output_files:
            print(f"Contigs already assembled for '{fasta_file}' with minimum overlap {min_overlap}")
            return assembled_contigs_file_path


    reads = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    print("Assembling reads...")
    print(f"Number of reads: {len(reads)}")
    print(f"Minimum overlap: {min_overlap}")

    # Find overlaps in parallel
    overlaps = find_overlaps_parallel(reads, min_overlap)

    # Build overlap graph
    graph = build_overlap_graph(overlaps)

    # Assemble genome in parallel
    contigs = assemble_genome_parallel(graph)
    print(f"Number of contigs before filtering: {len(contigs)}")
    filtered_contigs = filter_contigs(contigs)
    print(f"Number of contigs after filtering: {len(filtered_contigs)}")
    print(f"Longest contig length: {max(len(filtered_contigs) for filtered_contigs in filtered_contigs)}")

    with open(assembled_contigs_file_path, "w") as f:
        for i, contig in enumerate(filtered_contigs):
            f.write(f">contig_{i + 1}\n{contig}\n")

    print(f"Assembled contigs saved to '{assembled_contigs_file_path}'\n\n")
    return assembled_contigs_file_path
