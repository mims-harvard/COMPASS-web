#! /bin/bash

#sbatch --mem 100G -c 8 -t 5-12:00 -p priority  ./tcga_download_process.sh

Rscript 01_tcga_download_process.R
