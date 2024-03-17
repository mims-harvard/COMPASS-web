#! /bin/bash


#>>>02_fastqc_multiqc.sh>>>
work_dir=/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/fastaq/Gide
#fastqc
fastqc -t 20 -o ${work_dir}/data/rawdata/ ${work_dir}/data/rawdata/ERR*.fastq.gz
#multiqc
multiqc -o ${work_dir}/data/rawdata/ ${work_dir}/data/rawdata/*.zip

#>>>02_fastqc_multiqc.sh>>>


# sbatch --mem 32G -c 20 -t 5-12:00 -p priority ./02_fastqc_multiqc.sh
