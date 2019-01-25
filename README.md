# About

This repository contains the data analysis code used in the manuscript:

**Genetic effects on promoter usage are highly context-specific and contribute to complex traits**
https://doi.org/10.7554/eLife.41673

## Directories
* **analysis** - data analysis code
* **configs** - Snakemake workflow config files
* **scripts** - Stand-alone scripts used by the Snakemake workflows

## Snakemake workflows
* **quantify_transcription.snakefile** - RNA-seq quantification workflow
* **map_QTLs.snakefile** - QTL mapping workflow
* **run_coloc.snakefile** - Colocalisation workflow
* **run_fgwas.snakefile** - fgwas enrichment workflow
* **quantify_simulations.snakefile** - Quantify simulated RNA-seq data

## Dependencies
* **txrevise** - https://github.com/kauralasoo/txrevise
* Many functions are imported from the seqUtils package: https://github.com/kauralasoo/seqUtils
