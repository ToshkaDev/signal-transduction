# SIGNAL: A high-throughput pipeline for large-scale analysis of microbial signal transduction systems

![License](https://img.shields.io/github/license/ToshkaDev/signal-transduction)
![Python](https://img.shields.io/badge/python-3.6%2B-blue)
![Shell](https://img.shields.io/badge/shell-Bash-blue)
![Repo Size](https://img.shields.io/github/repo-size/ToshkaDev/signal-transduction)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)

> 🚧 **This project is under active development. Expect frequent changes.**

**SIGNAL** (**S**ystematic **I**nvestigation of **G**enomic **N**etworks for **A**nalysis of **L**ogic-based signaling) is a pipeline designed to perform high-throughput analysis of bacterial and archaeal genomes to uncover patterns in signal transduction systems across genomes, taxonomic groups, and functional architectures. The current dataset includes **26,221 representative genomes**.

---

## System Requirements

### Hardware

A standard computer with sufficient RAM for in-memory processing should be adequate.

### Software

- Python 3.6 or higher

### Tested Operating Systems

- Ubuntu 20.04 
- Linux Mint 20.2

---

## Installation

Clone the repository using Git:

``` git clone https://github.com/ToshkaDev/signal-transduction.git ```

This will download the repository and set up the pipeline for use.

### Usage

To launch the pipeline, use the master script:

```
cd signal-transduction
./analyze.sh
```

The script will first check whether the initial long-running step has already been completed by examining the presence of files in **results/obtain_and_process_tcs/**. Based on this, it will either start the entire pipeline or skip the completed steps.

## Pipeline Steps

### 1. Preparation

- Unpacks archived input files  
- Extracts genome lists (bacterial and archaeal)  
- Assigns genome sources (MiST genomes or MiST MAGs databases)  
- Creates all necessary input and output directories

Input files include:
- A dataset of 26,221 bacterial and archaeal genomes. It is also possible to use your own list of genomes prepared in accordance with the format used.
- Signal transduction domain definitions from the <a href="https://mistdb.com/" target="_blank">MiST (Microbial Signal Transduction)</a> database
- A metadata file from the <a href="https://gtdb.ecogenomic.org/" target="_blank">Genome Taxonomy Database (GTDB, release r214)</a>

### 2. TCS Extraction and Processing (`obtain_and_process_tcs.py`)

- Fetches two-component systems from the MiST database using its API  
- Analyzes protein domain compositions and architectures  
- Outputs tabulated results listing:
  - Genomes  
  - Histidine kinases (HKs) and response regulators (RRs)  
  - Their protein domain compositions and architectures

### 3. Genome-Level Analysis (`analyze_tcs_per_genome.py`)

- Analyzes and reports domain composition statistics for HKs and RRs per genome  
- Reports:
  - Number and type of input domains in HKs  
  - Additional domains in RRs  
- Normalizes statistics by genome size and total number of encoded proteins

### 4. Taxonomy-Level Analysis (`analyze_tcs_per_taxon.py`)

- Analyzes domain composition statistics for HKs and RRs at each taxonomic level:
  - Species  
  - Genus  
  - Family  
  - Order  
  - Class  
  - Phylum  
  - Kingdom  
- Normalizes results by the number of genomes per taxonomic level

The <a href="https://gtdb.ecogenomic.org/" target="_blank">GTDB</a> taxonomy is used.


## Future Work

- Visualization modules for domain architecture patterns


