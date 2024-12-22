import os
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner

def detect_misassemblies(contigs, reference, breakpoint_threshold=1000, mismatch_threshold=10):
    misassemblies = 0
    assembly_length = sum(len(contig) for contig in contigs)

    # Initialize PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2

    for contig in contigs:
        alignments = aligner.align(reference, contig)
        if not alignments:
            continue

        # Analyze all alignments
        for alignment in alignments:
            ref_coords, contig_coords = alignment.aligned
            for (ref_start, ref_end), (contig_start, contig_end) in zip(ref_coords, contig_coords):
                # Check for misassemblies
                gap_length = (ref_end - ref_start) - (contig_end - contig_start)
                mismatches = alignment.score / aligner.match_score
                if gap_length > breakpoint_threshold or mismatches > mismatch_threshold:
                    misassemblies += 1


    misassembly_rate = calculate_misassembly_rate(contigs, misassemblies, assembly_length)
    return misassemblies, misassembly_rate

def calculate_misassembly_rate(contigs, misassemblies, assembly_length):
    assembly_length_100kb = assembly_length / 100000
    misassembly_rate = misassemblies / assembly_length_100kb
    return misassembly_rate

def calculate_n50(contigs):
    total_length = sum(len(contig) for contig in contigs)
    sorted_contigs = sorted(contigs, key=len, reverse=True)
    cumulative_length = 0
    for contig in sorted_contigs:
        cumulative_length += len(contig)
        if cumulative_length >= total_length / 2:
            return len(contig)
    return 0


def calculate_genome_fraction(assembly, reference):
    assembly_seq = ''.join(assembly)
    covered_bases = 0
    for i in range(len(reference)):
        if i < len(assembly_seq) and reference[i] == assembly_seq[i]:
            covered_bases += 1
    return (covered_bases / len(reference)) * 100


def generate_assembly_report(error_free_file, error_prone_file, reference_genome_file):
    # Load contigs and reference genome
    error_free_contigs = list(SeqIO.parse(error_free_file, "fasta"))
    error_prone_contigs = list(SeqIO.parse(error_prone_file, "fasta"))
    reference_genome = str(next(SeqIO.parse(reference_genome_file, "fasta")).seq)

    # Compute metrics
    metrics = defaultdict(dict)
    for assembly_type, contigs in [("Error-free", error_free_contigs), ("Error-prone", error_prone_contigs)]:
        metrics[assembly_type]["Contig Count"] = len(contigs)
        metrics[assembly_type]["Largest Contig"] = max(len(contig) for contig in contigs)
        metrics[assembly_type]["N50"] = calculate_n50([str(contig.seq) for contig in contigs])
        metrics[assembly_type]["Genome Fraction"] = calculate_genome_fraction([str(contig.seq) for contig in contigs],
                                                                              reference_genome)
        metrics[assembly_type]["Total Assembly Length"] = sum(len(contig) for contig in contigs)
        metrics[assembly_type]["Misassemblies"], metrics[assembly_type]["Misassembly Rate"] = detect_misassemblies(contigs, reference_genome)

    # Generate report
    report = "Genome Assembly Report\n"
    report += "======================\n\n"

    report += "1. Assembly Statistics\n"
    report += "----------------------\n"
    for metric in ["Contig Count", "Largest Contig", "N50", "Genome Fraction", "Total Assembly Length", "Misassemblies", "Misassembly Rate"]:
        report += f"{metric}:\n"
        report += f"  Error-free: {metrics['Error-free'][metric]}\n"
        report += f"  Error-prone: {metrics['Error-prone'][metric]}\n\n"

    # Generate the histogram
    report += "2. Contig Length Distribution\n"
    report += "------------------------------\n"
    output_dir = "../data/output"

    # Prepare data for histograms
    error_free_lengths = [len(contig) for contig in error_free_contigs]
    error_prone_lengths = [len(contig) for contig in error_prone_contigs]

    # Plot the histograms
    plt.figure(figsize=(10, 6))
    plt.hist(error_free_lengths, bins=range(0, max(error_free_lengths + error_prone_lengths) + 100, 100),
             alpha=0.5, label="Error-free", density=True, color="green", edgecolor="black")
    plt.hist(error_prone_lengths, bins=range(0, max(error_free_lengths + error_prone_lengths) + 100, 100),
             alpha=0.5, label="Error-prone", density=True, color="red", edgecolor="black")

    # Add labels, legend, and title
    plt.xlabel("Contig Length")
    plt.ylabel("Relative Frequency")
    plt.title("Contig Length Distribution")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Save the figure
    output_path = os.path.join(output_dir, "contig_length_distribution.png")
    plt.savefig(output_path)
    plt.close()

    report += f"See '{output_path}' for a visual representation.\n\n"
    report += "3. Analysis\n"
    report += "------------\n"
    report += "Comparison between error-free and error-prone assemblies:\n"
    for metric in ["Contig Count", "Largest Contig", "N50", "Genome Fraction", "Total Assembly Length"]:
        diff = metrics["Error-free"][metric] - metrics["Error-prone"][metric]
        report += f"- {metric}: Difference of {diff}\n"

    report += "\nInterpretation:\n"
    report += "- A lower contig count generally indicates a more contiguous assembly.\n"
    report += "- Larger contigs and higher N50 values suggest better assembly quality.\n"
    report += "- Higher genome fraction indicates better coverage of the reference genome.\n"
    report += "- The total assembly length should ideally be close to the reference genome length.\n\n"

    report += "4. Conclusions\n"
    report += "--------------\n"
    if metrics["Error-free"]["N50"] > metrics["Error-prone"]["N50"]:
        report += "The error-free assembly appears to be more contiguous, as expected.\n"
    else:
        report += "Unexpectedly, the error-prone assembly shows higher contiguity. This might be due to aggressive merging of reads despite errors.\n"

    if metrics["Error-free"]["Genome Fraction"] > metrics["Error-prone"]["Genome Fraction"]:
        report += "The error-free assembly covers more of the reference genome, indicating better completeness.\n"
    else:
        report += "The error-prone assembly unexpectedly covers more of the reference genome. This might indicate over-assembly or duplication.\n"

    report += "\n5. Recommendations\n"
    report += "-------------------\n"
    report += "- Consider implementing error correction techniques to improve the error-prone assembly.\n"
    report += "- Experiment with different overlap lengths and mismatch tolerances in the assembly algorithm.\n"
    report += "- Validate assemblies using additional metrics such as BUSCO scores for gene completeness.\n"
    report += "- Consider using more sophisticated assembly algorithms that can handle repeats and resolve ambiguities better.\n"

    # Write report to file
    with open(f"{output_dir}\\assembly_report.txt", "w") as f:
        f.write(report)

    print("Report generated and saved as 'assembly_report.txt'")
    print("Contig length distribution plot saved as 'contig_length_distribution.png'")