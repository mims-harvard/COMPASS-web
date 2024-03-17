# RNA-seq Data Analysis Workflow

## 1. Quality Control with FastQC and MultiQC

- **FastQC**
  - Purpose: Checks the quality of raw sequencing data (typically FASTQ files).
  - Outputs: Quality reports detailing sequence quality scores, sequence length distribution, GC content, etc.

- **MultiQC**
  - Purpose: Aggregates results from multiple analyses (like FastQC) into a single report.
  - Outputs: A comprehensive quality control report for multiple samples.

## 2. Data Preprocessing with fastp

- **fastp**
  - Purpose: Cleans raw sequencing data to improve quality.
  - Tasks: Removes low-quality reads, trims adapters, filters low complexity sequences, etc.
  - Outputs: High-quality, cleaned sequencing data ready for alignment.

## 3. Alignment with STAR2 (Two-Pass Method)

- **STAR Aligner**
  - Purpose: Maps RNA-seq reads to a reference genome.
  - Special Feature: Two-pass method for enhanced detection of splice junctions.
  - Process:
    1. First Pass: Perform initial alignment.
    2. Second Pass: Utilize information from the first pass to better identify rare or novel splice junctions.
  - Outputs: Aligned sequencing reads to the reference genome, with improved accuracy for splice junctions.
