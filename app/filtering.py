import multiprocessing as mp
from Bio import SeqIO
from collections import defaultdict


def kmer_hash(sequence, k=10):
    return set(sequence[i:i + k] for i in range(len(sequence) - k + 1))


def is_subseq_optimized(smaller, larger, threshold, k=10):
    smaller_kmers = kmer_hash(smaller, k)
    larger_kmers = kmer_hash(larger, k)
    shared_kmers = len(smaller_kmers.intersection(larger_kmers))
    return shared_kmers / len(smaller_kmers) >= threshold


def filter_chunk(chunk, all_contigs, min_length, coverage_threshold):
    unique_contigs = []
    for contig in chunk:
        if len(contig.seq) < min_length:
            continue
        if not any(is_subseq_optimized(str(contig.seq), str(unique_contig.seq), coverage_threshold) for unique_contig in
                   unique_contigs):
            unique_contigs.append(contig)
    return unique_contigs


def filter_contigs(file_name, min_length=100, coverage_threshold=0.9):
    # Read contigs from file
    contigs = list(SeqIO.parse(file_name, "fasta"))

    # Sort contigs by length in descending order
    sorted_contigs = sorted(contigs, key=lambda x: len(x.seq), reverse=True)

    # Prepare chunks for parallel processing
    num_processes = mp.cpu_count()
    chunk_size = max(1, len(sorted_contigs) // num_processes)
    chunks = [sorted_contigs[i:i + chunk_size] for i in range(0, len(sorted_contigs), chunk_size)]

    # Process chunks in parallel
    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(filter_chunk,
                               [(chunk, sorted_contigs, min_length, coverage_threshold) for chunk in chunks])

    # Combine results
    filtered_contigs = [contig for chunk_result in results for contig in chunk_result]

    # Write filtered contigs back to the same file
    SeqIO.write(filtered_contigs, file_name, "fasta")

    print(f"Filtered contigs saved to {file_name}")
    print(f"Original contig count: {len(contigs)}")
    print(f"Filtered contig count: {len(filtered_contigs)}")

