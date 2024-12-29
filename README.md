# Genome Assembly Project

## Overview
This project focuses on genome assembly using synthetic reads generated from an existing reference genome. The synthetic reads are assembled under varying parameters, including different error probabilities, to investigate the effects of sequencing errors on assembly performance. By comparing assemblies generated from error-free reads and error-prone reads, this project aims to explore differences in contig length, assembly quality, genome fraction, and misassembly rates. 
The findings provide insights into how sequencing errors influence the accuracy and efficiency of genome assembly processes.

## Features
- Generate error-free and error-prone reads from a reference genome
- Assemble reads using an overlap graph approach
- Parallel processing for improved performance
- Contig filtering to remove redundant sequences
- Performance metrics calculation:
  - N50 statistic
  - Genome fraction
  - Misassembly rate
  - Contig count
  - Largest contig length
- Experimental parameter variation and analysis
- Visualization of assembly metrics


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
1. Configure parameters in config.yaml
   - Input Paths: Specify the locations of raw sequencing data.
   - Quality Thresholds: Set thresholds for trimming and filtering.
   - Assembly Parameters: Adjust settings for the assembler.
   - Output Directories: Define where results and logs should be stored.


2. Run the main experiment:
   ```bash
   python experimentation.py
   ```
   
3. For a single test run:
   ```bash
   python test.py
   ```

### Basic Experiment
```python
# Load your config parameters
config = load_config()

# Run experiments for all combinations of parameters
for n, l, p in itertools.product(n_values, l_values, p_values):
    result = run_experiment(n, l, p, genome_fasta_file, fasta_files_output_dir)
    
# Visualize findings
visualize_results(results_file, results_output_dir)
```

## Project Structure
```
Genome_Assembly/
├── app/
│   ├── generate_reads.py                         # Script to simulate or generate reads.
│   ├── assemble_genome.py                        # Main script for genome assembly.
│   ├── assemble_genome_improving_performance.py  # Performance-optimized genome assembly.
│   ├── filtering.py                              # Scripts to filter redundant contigs.
│   ├── generate_assembly_report.py               # Generates a detailed assembly report.
│
├── data/
│   ├── input/
│   │   ├── phiX_genome.fasta                     # Example input genome file.
│   │
│   ├── output/
│       ├── faste_files
│           ├── ...
│       ├── results
│           ├── contig_count_analysis.png
│           ├── experiment_results.txt
│           ├── largest_contig_analysis.png
│           ├── n50_analysis.png
│       ├── test_results
│           ├── assembled_contigs_l100_n10000_p0%.fasta   # Assembled contigs with specific parameters.
│           ├── assembled_contigs_l100_n10000_p1%.fasta   # Assembled contigs with different parameters.
│           ├── assembly_quality_performance_measures.pdf # PDF explanation of performance metrics.
│           ├── assembly_report.txt                     # Text-based assembly report.
│           ├── contig_length_distribution.png          # Visualization of contig length distribution.
│           ├── reads_l100_n10000_p0%.fasta             # Reads generated with specific parameters.
│           ├── reads_l100_n10000_p1%.fasta             # Reads generated with different parameters.
│
├── experimentation.py                        # For experimentation and testing of assembly methods.
├── test.py                                       # Entry point for running the pipeline.
├── config.yaml                                   # YAML configuration file for customizable parameters.
├── requirements.txt                              # Required Python dependencies for the project.
```


## Contact
Gilad Mevorach - [@giladmev](https://github.com/giladmev)

Project Link: [https://github.com/giladmev/Genome-Assembly](https://github.com/giladmev/Genome-Assembly)

