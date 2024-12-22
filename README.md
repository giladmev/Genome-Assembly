# Genome Assembly Project

## Overview
This project focuses on genome assembly using synthetic reads generated from an existing reference genome. The synthetic reads are assembled under varying parameters, including different error probabilities, to investigate the effects of sequencing errors on assembly performance. By comparing assemblies generated from error-free reads and error-prone reads, this project aims to explore differences in contig length, assembly quality, genome fraction, and misassembly rates. 
The findings provide insights into how sequencing errors influence the accuracy and efficiency of genome assembly processes.

## Features
- Overlap graph construction and genome assembly algorithm implementation
- Customizable parameters (k-mer length, error rate, coverage)
- Performance metrics calculation:
  - N50 statistic
  - Genome fraction
  - Misassembly rate
  - Contig count
  - Largest contig length
- Experimental parameter variation and analysis
- Visualization tools for assembly results

## Requirements
- Python 3.8+
- Required Python packages:
  ```
  numpy
  matplotlib
  pandas
  biopython
  pyyaml
  ```

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/giladmev/Genome-Assembly.git
   cd Genome-Assembly
   ```

2. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
### Configuration
The config.yaml file allows customization of various parameters:

- Input Paths: Specify the locations of raw sequencing data.
- Quality Thresholds: Set thresholds for trimming and filtering.
- Assembly Parameters: Adjust settings for the assembler.
- Output Directories: Define where results and logs should be stored.

### Basic Usage
```python
# Load your config parameters
config = load_config()

# Generate your sequencing reads
reads = generate_reads("path_to_reads.fasta")

# Run assembly
contigs = assemble_reads(reads_file_path)

# Calculate assembly metrics and generate report
generate_assembly_report(contigs_file_path)
```

## Project Structure
```
Genome_Assembly/
├── app/
│   ├── generate_reads.py                         # Script to simulate or generate reads.
│   ├── assemble_genome.py                        # Main script for genome assembly.
│   ├── assemble_genome_improving_performance.py  # Performance-optimized genome assembly.
│   ├── filtering.py                              # Scripts to filter redundant contigs.
│   ├── experimentation.py                        # For experimentation and testing of assembly methods.
│   ├── generate_assembly_report.py               # Generates a detailed assembly report.
│
├── data/
│   ├── input/
│   │   ├── phiX_genome.fasta                     # Example input genome file.
│   │
│   ├── output/
│       ├── assembled_contigs_l100_n10000_p0%.fasta   # Assembled contigs with specific parameters.
│       ├── assembled_contigs_l100_n10000_p1%.fasta   # Assembled contigs with different parameters.
│       ├── assembly_quality_performance_measures.pdf # PDF explanation of performance metrics.
│       ├── assembly_report.txt                     # Text-based assembly report.
│       ├── contig_length_distribution.png          # Visualization of contig length distribution.
│       ├── reads_l100_n10000_p0%.fasta             # Reads generated with specific parameters.
│       ├── reads_l100_n10000_p1%.fasta             # Reads generated with different parameters.
│
├── main.py                                       # Entry point for running the pipeline.
├── config.yaml                                   # YAML configuration file for customizable parameters.
├── requirements.txt                              # Required Python dependencies for the project.
```


## Contact
Gilad Mevorach - [@giladmev](https://github.com/giladmev)

Project Link: [https://github.com/giladmev/Genome-Assembly](https://github.com/giladmev/Genome-Assembly)

