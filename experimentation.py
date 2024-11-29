import itertools
from generate_reads import generate_reads
from assamble_genome import assemble_reads
from generate_assembly_report import calculate_n50

def run_experiment(n, l, p):
    # Generate reads with specified parameters
    generate_reads("data\\input\\phiX_genome.fasta", l, n, error_prob=p)

    # Assemble reads and compute metrics
    contigs = assemble_reads("data\\output\\error_prone_reads.fasta", min_overlap=20)
    n50 = calculate_n50(contigs)
    contig_count = len(contigs)
    largest_contig = max(len(contig) for contig in contigs)

    # Print results
    print(f"Parameters: n={n}, l={l}, p={p}")
    print(f"N50: {n50}, Contig Count: {contig_count}, Largest Contig: {largest_contig}\n")

    return {
        "n": n,
        "l": l,
        "p": p,
        "N50": n50,
        "Contig Count": contig_count,
        "Largest Contig": largest_contig
    }


# Define parameter ranges
n_values = [5000, 10000, 20000]
l_values = [50, 100, 150]
p_values = [0, 0.01, 0.05]

# Run experiments for all combinations of parameters
results = []
for n, l, p in itertools.product(n_values, l_values, p_values):
    result = run_experiment(n, l, p)
    results.append(result)

# Document findings
with open("data\\output\\experiment_results.txt", "w") as f:
    f.write("n\tl\tp\tN50\tContig Count\tLargest Contig\n")
    for result in results:
        f.write(
            f"{result['n']}\t{result['l']}\t{result['p']}\t{result['N50']}\t{result['Contig Count']}\t{result['Largest Contig']}\n")

print("Experiment results saved to 'experiment_results.txt'")