## Welcome to compass pages

## **compass 101**
#### This page provides details about the data pre-processing for the Compass input.


## COMPASS data preprocessing 

[![Codebase](https://img.shields.io/badge/Codebase-Github-green)](https://github.com/mims-harvard/compass)
[![Slack](https://img.shields.io/badge/Project-Slack-orange)](https://zitniklab-harvard.slack.com/archives/C05S6LEQ3ED)



Data Processing Pipeline Repository
This repository provides a pipeline for generating TPM (Transcripts Per Million) expression data from FASTQ files, as well as scripts and examples for processing TCGA (The Cancer Genome Atlas) data. The following sections outline the repository structure, usage instructions, and important considerations.

High-throughput sequencing (RNA-seq) data often require a series of preprocessing and normalization steps before downstream analyses. These steps typically include:

Quality control and filtering of raw reads from FASTQ files.
Alignment to a reference genome or transcriptome.
Quantification of gene or transcript expression (commonly expressed as TPM, RPKM, FPKM, etc.).
Processing of data from public databases like TCGA (The Cancer Genome Atlas), which can involve downloading, reformatting, and integrating with other data for comprehensive analysis.
This repository aims to provide a reproducible pipeline to guide users from raw FASTQ data to TPM expression results, as well as an example workflow to handle TCGA data.

