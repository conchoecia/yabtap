# yabtap - Yet another bulk transcriptome assembly pipeline

This snakemake pipeline is useful for managing transcriptome assemblies over many sequencing runs. 
It takes a config file defining the samples, the libraries that belong to those samples, and any associated metadata. 
Using the metadata it assembles the transcriptomes using Trinity (and eventually also RNAspades), translates the output, automatically formats and renames the sequences, and makes nucleotide blast, protein blast, and diamond databases. The pipeline also performs transcript quantification of each library using Kallisto.

## Benefits to using BATs

The benefits to assembling transcriptomes using BATs, as opposed to assembling them individually on a case-by-case basis:
  - Hands-on time is minimal. Fill out new entries in the config file, and all processing steps are completed automatically.
  - It is reproducible. Trimming and assembly parameters are consistent across samples. This is useful for large-scale comparative studies such as phylogenomics.
  - It is computationally efficient. All steps are handled by Snakemake to ensure that steps are only run when necessary, and are parallelized to finish the assemblies faster than a for-loop based pipeline.
  - It can be restarted easily. If a run crashes, simply execute the snakemake command again and the pipeline will pick up where it left off.
  - It has useful additional features. BATs provides transcript counts and makes databases, saving extra hassle.

## Requirements

The software required for this pipeline:
- python3 with biopython
- snakemake
- Docker
- kallisto
- bioawk
- Transdecoder

## Parameters

- Trimming
  - Adapters: The adapters used in the trimming protocol are a collection of all of the adapters found in the Trimmomatic software, and additional adapters for the Illumina smallRNA library prep kit. This pipeline aggressively trims adatpers and favors false positives over allowing some reads with adapters to pass into the normalization and de novo assembly steps.
