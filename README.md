# COMPASS: Generalizable AI predicts immunotherapy outcomes across cancers and treatments

[![ProjectPage](https://img.shields.io/badge/ProjectPage-COMPASSWebsite-red)](https://www.immuno-compass.com/)
[![Codebase](https://img.shields.io/badge/ProjectCode-COMPASSCode-green)](https://github.com/mims-harvard/COMPASS)
[![Slack](https://img.shields.io/badge/Project-Slack-orange)](https://zitniklab-harvard.slack.com/archives/C05S6LEQ3ED)

This repository provides details about the data pre-processing for the COMPASS input. It provides a reproducible pipeline to guide users from raw FASTQ data to TPM expression results, as well as an example workflow to handle TCGA data.

## COMPASS data processing pipeline repository

This repository provides a pipeline for generating TPM (Transcripts Per Million) expression data from FASTQ files, as well as scripts and examples for processing TCGA (The Cancer Genome Atlas) data. The following sections outline the repository structure, usage instructions, and important considerations. RNA-seq data often require a series of preprocessing and normalization steps before downstream analyses. These steps typically include:
* Quality control and filtering of raw reads from FASTQ files.
* Alignment to a reference genome or transcriptome.
* Quantification of gene or transcript expression (commonly expressed as TPM, RPKM, FPKM, etc.).
* Processing of data from public databases like TCGA (The Cancer Genome Atlas), which can involve downloading, reformatting, and integrating with other data for comprehensive analysis.


