# sanger_targseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.12.0-edge-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)](https://docs.docker.com/)

A comprehensive Nextflow pipeline for automated analysis of Sanger sequencing data, from raw .ab1 trace files to pathogen identification and phylogenetic analysis.

## 🔬 Overview

This pipeline streamlines the analysis of Sanger sequencing data for veterinary pathogen detection, providing:
- **Consensus assembly** from .ab1 trace files
- **Pathogen typing** using custom databases
- **BLAST-based identification** 
- **Phylogenetic analysis** with multiple sequence alignment and tree construction
- **Comprehensive HTML reports** with visualizations

## ✨ Features

- **📊 Tracy Assembly**: Assembles .ab1 trace files into high-quality consensus sequences
- **🔍 BLAST Analysis**: Performs searches on both consensus sequences and individual reads
- **🧬 Pathogen Typing**: Custom abricate database screening for selected veterinary pathogen detection and serotyping
- **🌳 Phylogenetic Analysis**: Multiple sequence alignment (MAFFT) and phylogenetic tree construction (IQ-TREE)
- **📈 Automated Reporting**: Generates comprehensive HTML reports with results and visualizations


## 🚀 Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥22.12.0-edge)
- [Docker](https://docs.docker.com/get-docker/)
- WSL2 (for Windows users)

### Basic Usage

```bash
# Clone the repository
git clone https://github.com/dhineshp565/sanger_targseq.git
cd sanger_targseq

# Run the pipeline
nextflow run main.nf --input /path/to/input --out_dir /path/to/output
```

### Input Structure

Your input directory should contain subdirectories, each with .ab1 trace files and reads files in txt or seq format:

```
input_directory/
├── sample1/
│   ├── trace1.ab1
│   ├── trace1.txt          
│   ├── trace2.ab1
│   ├── trace2.seq          
│   └── trace3.ab1
├── sample2/
│   ├── trace1.ab1
│   ├── trace1.txt          
│   ├── trace2.ab1
│   └── trace2.txt          
└── sample3/
    ├── trace1.ab1
    └── trace1.seq         
```

**File Types:**
- **Required**: `.ab1` - Sanger sequencing trace files for consensus assembly
                `.txt` or `.seq` - Associated reads files containing sequence data

## ⚙️ Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Path to input directory containing subdirectories with .ab1 files | `/data/sanger_data/` |
| `--out_dir` | Output directory path | `/results/sanger_analysis/` |
| `--blastdb_path` | Path to BLAST database | `/data/referenceDB/blast/microbe_db` |
| `--blastdb_name` | Name of BLAST database | `microbe_db` |


### Example Commands

```bash
# Basic run
nextflow run main.nf --input /data/samples --out_dir /results

# With custom BLAST database
nextflow run main.nf \
  --input /data/samples \
  --out_dir /results \
  --blastdb_path /custom/blast/db \
  --blastdb_name custom_db

```

## 📁 Output Structure

```
output_directory/
├── tracy_assembly/          # Consensus sequences and alignment files
├── abricate/               # Pathogen typing results
├── blast_consensus/        # BLAST results for consensus sequences
├── blast_sanger_reads/     # BLAST results for individual reads
├── orfipy/                # ORF predictions
├── seq_length/            # Sequence length statistics
├── mafft/                 # Multiple sequence alignments
├── iqtree/                # Phylogenetic trees
├── ggtree/                # Tree visualizations
├── LIMS/                  # Laboratory information management files
├── execution/             # Pipeline execution reports
└── sanger_targseq_report.html  # Main results report
```

## 🛠️ Pipeline Workflow

1. **Input Processing**: Creates sample list from input directory structure
2. **Tracy Assembly**: Assembles .ab1 files into consensus sequences
3. **Quality Control**: Generates reverse complement sequences for analysis
4. **Pathogen Typing**: Screens sequences against veterinary pathogen database
5. **BLAST Analysis**: Searches against reference databases
6. **ORF Prediction**: Identifies open reading frames
7. **Phylogenetic Analysis**: Multiple alignment and tree construction
8. **Reporting**: Generates comprehensive HTML reports with visualizations

## 🔧 Dependencies

### Core Tools (Containerized)

- **[Tracy](https://github.com/gear-genomics/tracy)** (v0.7.8): Sanger trace file assembly
- **[Abricate](https://github.com/tseemann/abricate)** (v1.0.1): Antimicrobial resistance gene detection
- **[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)** (v2.16.0): Sequence similarity search
- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** (v7.526): Multiple sequence alignment
- **[IQ-TREE](http://www.iqtree.org/)** (v2.4.0): Phylogenetic tree reconstruction
- **[OrfIpy](https://github.com/urmi-21/orfipy)** (v0.0.4): ORF prediction
- **[R/RMarkdown](https://rmarkdown.rstudio.com/)** (v2.10): Report generation
- **[seqtk](https://github.com/lh3/seqtk)** (v1.4): Sequence processing utilities


## 🐛 Troubleshooting

### Common Issues

1. **Input directory structure**: Verify subdirectories contain .ab1 files
2. **Permission errors**: Check file/directory permissions


## �️ Software Used

This pipeline integrates the following software tools:

### Core Analysis Tools
- **[Tracy](https://github.com/gear-genomics/tracy)** (v0.7.8) - Sanger trace file basecalling, alignment, and assembly
- **[Abricate](https://github.com/tseemann/abricate)** (v1.0.1) - Mass screening of contigs for antimicrobial resistance and virulence genes
- **[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)** (v2.16.0) - Basic Local Alignment Search Tool for sequence similarity searches
- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** (v7.526) - Multiple sequence alignment program
- **[IQ-TREE](http://www.iqtree.org/)** (v2.4.0) - Efficient phylogenomic software for maximum likelihood phylogenetic inference
- **[seqtk](https://github.com/lh3/seqtk)** (v1.4) - Toolkit for processing sequences in FASTA/FASTQ formats

### Additional Tools
- **[OrfIpy](https://github.com/urmi-21/orfipy)** (v0.0.4) - Fast and flexible tool for extracting ORFs from nucleotide sequences
- **[R](https://www.r-project.org/)** with **[RMarkdown](https://rmarkdown.rstudio.com/)** (v2.10) - Statistical computing and dynamic report generation
- **[ggtree](https://guangchuangyu.github.io/software/ggtree/)** - Phylogenetic tree visualization and annotation

### Workflow Management
- **[Nextflow](https://www.nextflow.io/)** (≥22.12.0-edge) - Workflow management system
- **[Docker](https://www.docker.com/)** - Containerization platform for reproducible environments

### Container Images
All tools are containerized using Docker images from:
- Quay.io Biocontainers
- Docker Hub official repositories  
- StaPH-B (State Public Health Bioinformatics) collection
- NCBI official images

*Note: Exact container versions and sources can be found in the `config/docker.config` file.*

## 📝 License

This project is open source and available under the [MIT License](LICENSE).

## 👤 Author

**dhineshp565**
- GitHub: [@dhineshp565](https://github.com/dhineshp565)
- Pipeline: [sanger_targseq](https://github.com/dhineshp565/sanger_targseq)

---

**Version**: v1.0.6  
**Last Updated**: September 2025
