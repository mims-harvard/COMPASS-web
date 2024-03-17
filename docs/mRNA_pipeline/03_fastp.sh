#! /bin/bash


#>>>fastp.sh>>>
work_dir=/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/fastaq/Gide
rawdata_dir=${work_dir}/data/rawdata
cleandata_dir=${work_dir}/data/cleandata/fastp
n_jobs=40
#cut -f 4 ${rawdata_dir}/filereport_read_run_PRJNA229998_tsv.txt | sed '1d' > ${rawdata_dir}/run_accession.txt

cd ${rawdata_dir}
ls *.fastq.gz | awk -F'[_.]' '{print $1}' | sort -u > run_accession.txt

## split into 3 files, sbatch each task
split -n l/5 -d ${rawdata_dir}/run_accession.txt ${rawdata_dir}/run_accession_part

cat ${rawdata_dir}/run_accession.txt | while read i
do
fastp \
--in1 ${rawdata_dir}/${i}_1.fastq.gz \
--in2 ${rawdata_dir}/${i}_2.fastq.gz  \
--out1 ${cleandata_dir}/${i}_1.fastp.fq.gz \
--out2 ${cleandata_dir}/${i}_2.fastp.fq.gz \
--json ${cleandata_dir}/${i}.fastp.json \
--html ${cleandata_dir}/${i}.fastp.html \
--report_title ${cleandata_dir}/${i} \
--thread ${n_jobs}
done
#<<<fastp.sh<<<



# sbatch --mem 24G -c 40 -t 5-12:00 -p priority ./fastp.sh
