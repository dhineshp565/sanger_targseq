# sanger_targseq
This pipeline automates the analysis of Sanger sequencing data, streamlining the process from raw .ab1 trace files to BLAST search.


## Features

* Assembles .ab1 trace files into consensus sequences using tracy.
* Performs BLAST searches on:
    -Assembled consensus sequences
    -Forward and reverse individual reads
* Pathogen Typing with Abricate - Screens consensus sequences against a custom database for veterinary pathogen detection and serotyping.
* Taxonomic Classification using kraken2 - Assigns taxonomic information to consensus sequences


Requires input directory with sub-directories with ab1 trace files
Outputs a html report  with abricate and BLAST result, directoreies with resulst from each processes


Usage:
```
nextflow run main.nf --input path_to_input --out_dir Path_to_output
```

```
Parameters:

--input      Path to input directory
--out_dir    Output directory

Optional
--kraken_db     Path to a kraken database
--blastdb_path  Path to a BLAST database
--blastdb_name  Name of the BLAST database
	
```
## Dependencies
* nextflow
* docker
* wsl2
## Software Used - Docker image versions can be found in the docker.config file
* Tracy (Rausch, T., Fritz, M.H., Untergasser, A. and Benes, V. Tracy: basecalling, alignment, assembly and deconvolution of sanger chromatogram trace files. BMC Genomics 21, 230 (2020).https://doi.org/10.1186/s12864-020-6635-8)
* Abricate (https://github.com/tseemann/abricate)
* BLAST+ (Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990 Oct 5;215(3):403-10. doi: 10.1016/S0022-2836(05)80360-2. PMID: 2231712. https://github.com/ncbi/blast_plus_docs)
* rmarkdown (https://rmarkdown.rstudio.com/)
