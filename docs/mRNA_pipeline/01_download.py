# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 16:40:27 2023

@author: Wanxiang Shen
"""



#nohup python -u 01_download.py > download.out 2>&1 &
#download from https://www.ebi.ac.uk/ena/

import os
from tqdm import tqdm
import pandas as pd
from joblib import delayed,Parallel
from itertools import chain


cohort_folder = '/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/fastaq/Gide/'
rawdata_folder = os.path.join(cohort_folder, 'data', 'rawdata')
fastp_folder = os.path.join(cohort_folder, 'data', 'cleandata', 'fastp')
star_folder = os.path.join(cohort_folder, 'STAR')

os.system('mkdir -p %s' % rawdata_folder)
os.system('mkdir -p %s' % fastp_folder)
os.system('mkdir -p %s' % star_folder)

url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB23709&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_md5,fastq_ftp,submitted_md5,submitted_ftp,sra_md5,sra_ftp,sample_title&format=tsv&download=true&limit=0'
df = pd.read_csv(url, sep='\t')
df.to_csv(os.path.join(rawdata_folder,'PRJEB23709.csv'))

md5list = df.fastq_md5.apply(lambda x:x.split(';')).tolist()
md5list = list(chain(*md5list))

urlist = df.fastq_ftp.apply(lambda x:x.replace('ftp.sra.ebi.ac.uk/','fasp.sra.ebi.ac.uk:')).apply(lambda x:x.split(';')).to_list()
urlist = list(chain(*urlist))
dfd = pd.DataFrame([md5list, urlist]).T
dfd.columns=['md5', 'url']
dfd['dist_filename'] = dfd.url.apply(os.path.basename).apply(lambda x:os.path.join(rawdata_folder,x))

cmdf = lambda x: "ascp -i ~/micromamba/envs/RNA/etc/asperaweb_id_dsa.openssh -QT -l 300m -P33001 era-fasp@%s %s" % (x.url,x.dist_filename)
dfd['cmds'] = dfd.apply(cmdf, axis=1)

def download(cmd):
    res = os.system(cmd)
    if res ==0:
        print('%s: Success!' % cmd)     
    else:
        print('%s: Failed!' % cmd)
    return res


cmds = dfd.cmds.to_list()
P = Parallel(n_jobs=10)
res = P(delayed(download)(cmd) for cmd in tqdm(cmds, ascii=True)) 



dfd.to_csv(os.path.join(rawdata_folder, 'file.meta.csv'))

print('****************DONE**********************')