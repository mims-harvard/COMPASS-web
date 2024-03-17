import pandas as pd
from glob import glob
import os

from rnanorm.datasets import load_toy_data
from rnanorm import FPKM, TPM, CPM, TMM 

# load the gene length
gtf_path = '/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/GENCODE/release_36/HS.gencode.v36.annotation.gtf'
tpm = TPM(gtf_path).set_output(transform="pandas")
tmm = TMM(m_trim = 0.3, a_trim = 0.05).set_output(transform="pandas")


out_dir='/n/data1/hms/dbmi/zitnik/lab/users/was966/ITRP/fastaq/Gide/STAR/'
outs = glob(os.path.join(out_dir, "*ReadsPerGene.out.tab"))


counts=[]
for out in outs:
    ids = os.path.basename(out).replace('ReadsPerGene.out.tab', '')
    df = pd.read_csv(out,sep='\t', header=None)
    df.columns=['geneid', 'unstranded_counts', 'unstranded_counts_1st','unstranded_counts_2nd']
    df = df.iloc[4:].set_index('geneid')
    cts = df.unstranded_counts
    cts.name = ids
    counts.append(cts)

## counts, TPM, TMM
df_counts = pd.concat(counts,axis=1).T
df_tpm = tpm.fit_transform(df_counts)
df_tmm = tmm.fit_transform(df_counts)

df_counts.to_csv('./raw/untranded_counts.csv')
df_tpm.to_csv('./raw/unstranded_tpm.csv')