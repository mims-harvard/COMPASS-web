#!/bin/sh

#>>>star.sh>>>

work_dir=/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/fastaq/Gide
index_dir=/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/star_index/gdc_index/star-2.7.5c_GRCh38.d1.vd1_gencode.v36

rawdata_dir=${work_dir}/data/rawdata
input_dir=${work_dir}/data/cleandata/fastp
out_dir=${work_dir}/STAR
n_jobs=20


cat ${rawdata_dir}/run_accession_part04 | while read i
do
STAR \
--readFilesIn ${input_dir}/${i}_1.fastp.fq.gz ${input_dir}/${i}_2.fastp.fq.gz \
--outSAMattrRGline ID:sample SM:sample PL:ILLUMINA \
--genomeDir ${index_dir} \
--readFilesCommand zcat \
--runThreadN ${n_jobs} \
--twopassMode Basic \
--outFileNamePrefix ${out_dir}/${i} \
--quantMode GeneCounts \
--outSAMtype BAM Unsorted \
--outSAMunmapped None \

done


#<<<star.sh<<<
#sbatch --mem 64G -c 20 -t 29-12:00 -p long ./04.star_align_counts_04.sh
