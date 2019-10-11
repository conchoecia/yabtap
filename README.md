# BATs - Bulk Assembly of Transcriptomes

This snakemake pipeline is useful for managing transcriptome assemblies over many sequencing runs. 
It takes a config file defining the samples, the libraries that belong to those samples, and any associated metadata. 
Using the metadata it assembles the transcriptomes using Trinity (and eventually also RNAspades), translates the output, automatically formats and renames the sequences, and makes nucleotide blast, protein blast, and diamond databases. The pipeline also performs transcript quantification of each library using Kallisto.

# Requirements

The software required for this pipeline:
- python3 with biopython
- snakemake
- Docker
- kallisto
- bioawk
- Transdecoder


