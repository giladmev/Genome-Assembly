# Genome Assembly Project

## Overview
This project implements a genome assembly pipeline using overlap graphs. It includes modules for generating error-free and error-prone reads, assembling reads into contigs, and analyzing assembly performance.

## Features
- Generate error-free and error-prone reads from a reference genome
- Assemble reads using an overlap graph approach
- Parallel processing for improved performance
- Contig filtering to remove redundant sequences
- Experimental parameter variation and analysis
- Visualization of assembly metrics


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

2. Run the main experiment:
   ```bash
   python experimentation.py
   ```
   
3. For a single test run:
   ```bash
   python test.py
   ```

### Modules
- generate_reads.py: Generate error-free and error-prone reads
- assamble_genome_improving_performance.py: Assemble reads into contigs
- filtering.py: Filter and remove redundant contigs
- experimentation.py: Run experiments with varying parameters
- test.py: Perform a single test run

### Configuration
Edit config.yaml to set:
- Input genome file
- Output directories
- Sequencing parameters (read length, number of reads)
- Assembly parameters (minimum overlap, error probability)

### Experimentation
The experimentation.py script:
- Runs assembly with different combinations of parameters
- Generates visualizations for N50, contig count, and largest contig size
- Saves results in TSV format and creates plots

### Performance Measures
- N50: Contig length at which 50% of the assembly is in contigs of this size or larger
- Contig Count: Total number of assembled contigs
- Largest Contig: Length of the longest assembled contig

### Results
Experiment results are saved in:
- experiment_results.txt: TSV file with metrics for each parameter combination
- PNG files in the results directory: Visualizations of metrics vs. parameters

## Contact
Gilad Mevorach - [@giladmev](https://github.com/giladmev)

Project Link: [https://github.com/giladmev/Genome-Assembly](https://github.com/giladmev/Genome-Assembly)

