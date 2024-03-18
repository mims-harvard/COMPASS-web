# Conceptor Input Requirement

Please note that the input for Conceptor should be mRNA's TPM (**Transcripts Per Million**) expression values, not raw counts or other forms of mRNA expression values. TPM calculation is similar to FPKM but differs in the normalization process. In TPM, all transcripts are normalized for length first. Then, instead of using the total overall read count for size normalization, the sum of the length-normalized transcript values is used as a size indicator.

Please fell free to contact me if you have any questions on this.

## Data Processing Recommendation

If your data is in raw sequence format (FASTQ) or as raw counts, we recommend processing it using the following bioinformatics pipeline. This recommendation is based on the fact that our pretrained TCGA (The Cancer Genome Atlas) data was processed using this pipeline, and using the same pipeline for your input data may yield better results.



### 1. [mRNA-Seq Alignment Workflow](https://github.com/mims-harvard/conceptor-101/tree/main/mRNA_pipeline)

The RNA-Seq Alignment Workflow follows these steps:

`fastqc/multiqc  --> fastp ---> STAR2 align (two-pass method)`

For more information, please refer to [GDC mRNA expression pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/).


![pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/images/RNA-Seq-DR32_Image.png)

#### Specific Process

1. Begin with quality control using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiqc](http://multiqc.info).
2. Proceed with data preprocessing using [fastp](https://github.com/OpenGene/fastp) to clean raw sequencing data and improve quality.
3. Finally, align the RNA-seq reads to a reference genome using STAR version 2.7.5c, which maps RNA-seq reads to the [reference genome](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files). While custom index files can be created, we use the reference genome files downloaded from GDC. The link for the specific reference genome file (star-2.7.5c_GRCh38.d1.vd1_gencode.v36.tgz) is available [here](https://api.gdc.cancer.gov/data/c0008693-0583-4eac-bd5c-583070763893).



### 2. [Converting mRNA Raw Counts to TPM](https://github.com/mims-harvard/conceptor-101/blob/main/mRNA_pipeline/05_counts2tpm.py)

If your data is in mRNA expression counts, you can convert the mRNA raw counts to TPM values using the following method. 
This process involves normalization using gene lengths, so you will need to download the gene annotation file (v36).on.gtf.gz


#### Step 1: Download the human GENCODE annotation file (v36)

Download the GENCODE human annotation file (version 36) from the following link:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
mv gencode.v36.annotation.gtf.gz ./data
```

#### step2: using rnanorm tool to convert Count to TPM

```bash
#install rnanorm
pip install rnanorm 
```
----
```python
from rnanorm import FPKM, TPM, CPM, TMM 
gtf_path = "./data/gencode.v36.annotation.gtf"
tpm = TPM(gtf_path).set_output(transform="pandas")
df_tpm = tpm.fit_transform(df_counts)

```

#### step3: Now let's test on an example file


```python
# import packages
import pandas as pd
from rnanorm import FPKM, TPM, CPM, TMM 
```


```python
# convert count to TPM based the gtf file
gtf_path = "./data/gencode.v36.annotation.gtf.gz"
tpm = TPM(gtf_path).set_output(transform="pandas")
```


```python
# example of the raw counts
df_counts = pd.read_csv('./data/toy_raw_counts.csv', index_col=0)
df_counts.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ENSG00000223972.5</th>
      <th>ENSG00000227232.5</th>
      <th>ENSG00000278267.1</th>
      <th>ENSG00000243485.5</th>
      <th>ENSG00000284332.1</th>
      <th>ENSG00000237613.2</th>
      <th>ENSG00000268020.3</th>
      <th>ENSG00000240361.2</th>
      <th>ENSG00000186092.6</th>
      <th>ENSG00000238009.6</th>
      <th>...</th>
      <th>ENSG00000198886.2</th>
      <th>ENSG00000210176.1</th>
      <th>ENSG00000210184.1</th>
      <th>ENSG00000210191.1</th>
      <th>ENSG00000198786.2</th>
      <th>ENSG00000198695.2</th>
      <th>ENSG00000210194.1</th>
      <th>ENSG00000198727.2</th>
      <th>ENSG00000210195.2</th>
      <th>ENSG00000210196.2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ERR2208944</th>
      <td>6</td>
      <td>201</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>1376</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>947</td>
      <td>178</td>
      <td>0</td>
      <td>582</td>
      <td>0</td>
      <td>4</td>
    </tr>
    <tr>
      <th>ERR2208928</th>
      <td>0</td>
      <td>222</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>2263</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2549</td>
      <td>459</td>
      <td>0</td>
      <td>1486</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ERR2208949</th>
      <td>1</td>
      <td>487</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>2544</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1783</td>
      <td>377</td>
      <td>0</td>
      <td>745</td>
      <td>0</td>
      <td>4</td>
    </tr>
    <tr>
      <th>ERR2208900</th>
      <td>13</td>
      <td>569</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>14</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>...</td>
      <td>13168</td>
      <td>4</td>
      <td>3</td>
      <td>1</td>
      <td>10988</td>
      <td>2702</td>
      <td>0</td>
      <td>3746</td>
      <td>0</td>
      <td>101</td>
    </tr>
    <tr>
      <th>ERR2208922</th>
      <td>0</td>
      <td>29</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>14029</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>5480</td>
      <td>1302</td>
      <td>0</td>
      <td>4160</td>
      <td>0</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 60660 columns</p>
</div>




```python
# example of the TPM values
df_tpm = tpm.fit_transform(df_counts)
df_tpm.to_csv('./data/toy_tpm.csv')
df_tpm.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ENSG00000223972.5</th>
      <th>ENSG00000227232.5</th>
      <th>ENSG00000278267.1</th>
      <th>ENSG00000243485.5</th>
      <th>ENSG00000284332.1</th>
      <th>ENSG00000237613.2</th>
      <th>ENSG00000268020.3</th>
      <th>ENSG00000240361.2</th>
      <th>ENSG00000186092.6</th>
      <th>ENSG00000238009.6</th>
      <th>...</th>
      <th>ENSG00000198886.2</th>
      <th>ENSG00000210176.1</th>
      <th>ENSG00000210184.1</th>
      <th>ENSG00000210191.1</th>
      <th>ENSG00000198786.2</th>
      <th>ENSG00000198695.2</th>
      <th>ENSG00000210194.1</th>
      <th>ENSG00000198727.2</th>
      <th>ENSG00000210195.2</th>
      <th>ENSG00000210196.2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ERR2208944</th>
      <td>0.460737</td>
      <td>19.821755</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.109294</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>...</td>
      <td>133.036440</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>69.629485</td>
      <td>45.171249</td>
      <td>0.0</td>
      <td>67.957710</td>
      <td>0.0</td>
      <td>7.837047</td>
    </tr>
    <tr>
      <th>ERR2208928</th>
      <td>0.000000</td>
      <td>25.382838</td>
      <td>2.271609</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>...</td>
      <td>253.675132</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>217.297235</td>
      <td>135.050420</td>
      <td>0.0</td>
      <td>201.175794</td>
      <td>0.0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>ERR2208949</th>
      <td>0.069637</td>
      <td>43.552412</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.297342</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>...</td>
      <td>223.052187</td>
      <td>1.751014</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>118.886282</td>
      <td>86.760220</td>
      <td>0.0</td>
      <td>78.887687</td>
      <td>0.0</td>
      <td>7.107055</td>
    </tr>
    <tr>
      <th>ERR2208900</th>
      <td>0.172087</td>
      <td>9.672981</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.263771</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.024656</td>
      <td>...</td>
      <td>219.469413</td>
      <td>1.331418</td>
      <td>1.167811</td>
      <td>0.323478</td>
      <td>139.272015</td>
      <td>118.203257</td>
      <td>0.0</td>
      <td>75.402463</td>
      <td>0.0</td>
      <td>34.112682</td>
    </tr>
    <tr>
      <th>ERR2208922</th>
      <td>0.000000</td>
      <td>2.360453</td>
      <td>1.617126</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>...</td>
      <td>1119.515705</td>
      <td>3.187378</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>332.563864</td>
      <td>272.712079</td>
      <td>0.0</td>
      <td>400.922453</td>
      <td>0.0</td>
      <td>9.702754</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 60660 columns</p>
</div>



### 3. Preparing the inputs for the conceptor

The Inputs of conceptor model including the cancer type information and TPM values, the genes are identified by gene name, and gene name can be mapped from a dictionary contains the gene Ensembl ID, Entrez gene ID, and gene name.

Please find the cancer code of your data from this table: [TCGA Study Abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations):
| Study Abbreviation | Study Name                                                       |
| ------------------ | ---------------------------------------------------------------- |
| LAML               | Acute Myeloid Leukemia                                           |
| ACC                | Adrenocortical carcinoma                                         |
| BLCA               | Bladder Urothelial Carcinoma                                     |
| LGG                | Brain Lower Grade Glioma                                         |
| BRCA               | Breast invasive carcinoma                                        |
| CESC               | Cervical squamous cell carcinoma and endocervical adenocarcinoma |
| CHOL               | Cholangiocarcinoma                                               |
| LCML               | Chronic Myelogenous Leukemia                                     |
| COAD               | Colon adenocarcinoma                                             |
| CNTL               | Controls                                                         |
| ESCA               | Esophageal carcinoma                                             |
| FPPP               | FFPE Pilot Phase II                                              |
| GBM                | Glioblastoma multiforme                                          |
| HNSC               | Head and Neck squamous cell carcinoma                            |
| KICH               | Kidney Chromophobe                                               |
| KIRC               | Kidney renal clear cell carcinoma                                |
| KIRP               | Kidney renal papillary cell carcinoma                            |
| LIHC               | Liver hepatocellular carcinoma                                   |
| LUAD               | Lung adenocarcinoma                                              |
| LUSC               | Lung squamous cell carcinoma                                     |
| DLBC               | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma                  |
| MESO               | Mesothelioma                                                     |
| MISC               | Miscellaneous                                                    |
| OV                 | Ovarian serous cystadenocarcinoma                                |
| PAAD               | Pancreatic adenocarcinoma                                        |
| PCPG               | Pheochromocytoma and Paraganglioma                               |
| PRAD               | Prostate adenocarcinoma                                          |
| READ               | Rectum adenocarcinoma                                            |
| SARC               | Sarcoma                                                          |
| SKCM               | Skin Cutaneous Melanoma                                          |
| STAD               | Stomach adenocarcinoma                                           |
| TGCT               | Testicular Germ Cell Tumors                                      |
| THYM               | Thymoma                                                          |
| THCA               | Thyroid carcinoma                                                |
| UCS                | Uterine Carcinosarcoma                                           |
| UCEC               | Uterine Corpus Endometrial Carcinoma                             |
| UVM                | Uveal Melanoma                                                   |

#### Step1. Add the cancer type information
Suppose your data are all from Melonoma, here is an example to generate the conceptor's cancer type


```python
df_cancer_type = pd.DataFrame([], index = df_counts.index)
df_cancer_type['cancer_type'] = 'SKCM'
df_cancer_type.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cancer_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ERR2208944</th>
      <td>SKCM</td>
    </tr>
    <tr>
      <th>ERR2208928</th>
      <td>SKCM</td>
    </tr>
    <tr>
      <th>ERR2208949</th>
      <td>SKCM</td>
    </tr>
    <tr>
      <th>ERR2208900</th>
      <td>SKCM</td>
    </tr>
    <tr>
      <th>ERR2208922</th>
      <td>SKCM</td>
    </tr>
  </tbody>
</table>
</div>



After that, we need to map the cancer type to cancer code:


```python
import json
with open('./data/cancer_code.json') as f:
    cancer_code_map = json.load(f)
df_cancer_type['cancer_type'] = df_cancer_type['cancer_type'].map(cancer_code_map)
df_cancer_type.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cancer_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ERR2208944</th>
      <td>25</td>
    </tr>
    <tr>
      <th>ERR2208928</th>
      <td>25</td>
    </tr>
    <tr>
      <th>ERR2208949</th>
      <td>25</td>
    </tr>
    <tr>
      <th>ERR2208900</th>
      <td>25</td>
    </tr>
    <tr>
      <th>ERR2208922</th>
      <td>25</td>
    </tr>
  </tbody>
</table>
</div>



#### Step2. Now lets map the df_counts to conceptor's input genes.
The dictionary below contains the gene Ensembl ID, Entrez gene ID, and gene name


```python
gene_map = pd.read_csv('./data/conceptor_gene_map.csv')
gene_map.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ensid</th>
      <th>gene_name</th>
      <th>ensid_v36</th>
      <th>gene_type</th>
      <th>gene_supertype</th>
      <th>entrezgene</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000121410</td>
      <td>A1BG</td>
      <td>ENSG00000121410.12</td>
      <td>protein_coding</td>
      <td>protein_coding</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000148584</td>
      <td>A1CF</td>
      <td>ENSG00000148584.15</td>
      <td>protein_coding</td>
      <td>protein_coding</td>
      <td>29974.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000175899</td>
      <td>A2M</td>
      <td>ENSG00000175899.15</td>
      <td>protein_coding</td>
      <td>protein_coding</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000166535</td>
      <td>A2ML1</td>
      <td>ENSG00000166535.20</td>
      <td>protein_coding</td>
      <td>protein_coding</td>
      <td>144568.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000128274</td>
      <td>A4GALT</td>
      <td>ENSG00000128274.17</td>
      <td>protein_coding</td>
      <td>protein_coding</td>
      <td>53947.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
df_tpm_conceptor = df_tpm[gene_map.ensid_v36]
df_tpm_conceptor.columns = df_tpm_conceptor.columns.map(gene_map.set_index('ensid_v36').gene_name)
df_tpm_conceptor.shape
```




    (25, 15672)




```python
df_tpm_conceptor.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A1BG</th>
      <th>A1CF</th>
      <th>A2M</th>
      <th>A2ML1</th>
      <th>A4GALT</th>
      <th>A4GNT</th>
      <th>AAAS</th>
      <th>AACS</th>
      <th>AADAC</th>
      <th>AADAT</th>
      <th>...</th>
      <th>ZWILCH</th>
      <th>ZWINT</th>
      <th>ZXDA</th>
      <th>ZXDB</th>
      <th>ZXDC</th>
      <th>ZYG11A</th>
      <th>ZYG11B</th>
      <th>ZYX</th>
      <th>ZZEF1</th>
      <th>ZZZ3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ERR2208944</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>859.203620</td>
      <td>73.019466</td>
      <td>11.942279</td>
      <td>1.947147</td>
      <td>86.527503</td>
      <td>9.236956</td>
      <td>3.918524</td>
      <td>17.763974</td>
      <td>...</td>
      <td>16.725820</td>
      <td>14.012390</td>
      <td>7.255890</td>
      <td>6.984732</td>
      <td>16.469669</td>
      <td>0.879873</td>
      <td>21.106355</td>
      <td>89.920944</td>
      <td>47.520979</td>
      <td>21.534480</td>
    </tr>
    <tr>
      <th>ERR2208928</th>
      <td>0.038627</td>
      <td>0.032171</td>
      <td>881.830260</td>
      <td>7.533515</td>
      <td>12.650118</td>
      <td>2.778540</td>
      <td>95.158662</td>
      <td>10.227978</td>
      <td>1.798357</td>
      <td>10.818061</td>
      <td>...</td>
      <td>34.613126</td>
      <td>42.500215</td>
      <td>12.806729</td>
      <td>12.317719</td>
      <td>18.357604</td>
      <td>0.526526</td>
      <td>35.163162</td>
      <td>60.709750</td>
      <td>52.413439</td>
      <td>25.859226</td>
    </tr>
    <tr>
      <th>ERR2208949</th>
      <td>0.030213</td>
      <td>0.012581</td>
      <td>504.984491</td>
      <td>50.836895</td>
      <td>5.900676</td>
      <td>0.611231</td>
      <td>106.174319</td>
      <td>8.090318</td>
      <td>4.960132</td>
      <td>28.557439</td>
      <td>...</td>
      <td>12.677251</td>
      <td>19.670726</td>
      <td>12.836934</td>
      <td>9.511444</td>
      <td>9.528438</td>
      <td>0.154435</td>
      <td>23.424874</td>
      <td>69.710920</td>
      <td>40.638326</td>
      <td>25.926391</td>
    </tr>
    <tr>
      <th>ERR2208900</th>
      <td>0.143579</td>
      <td>0.023916</td>
      <td>1940.416805</td>
      <td>0.182940</td>
      <td>4.014771</td>
      <td>0.813332</td>
      <td>46.225429</td>
      <td>2.235265</td>
      <td>0.042219</td>
      <td>40.342963</td>
      <td>...</td>
      <td>24.998045</td>
      <td>21.527292</td>
      <td>8.889713</td>
      <td>7.719297</td>
      <td>12.813737</td>
      <td>0.670318</td>
      <td>18.452489</td>
      <td>56.563242</td>
      <td>36.542147</td>
      <td>32.326600</td>
    </tr>
    <tr>
      <th>ERR2208922</th>
      <td>0.109992</td>
      <td>0.286277</td>
      <td>1534.682495</td>
      <td>0.860539</td>
      <td>11.887174</td>
      <td>0.061813</td>
      <td>83.537696</td>
      <td>5.875654</td>
      <td>1.482365</td>
      <td>16.994521</td>
      <td>...</td>
      <td>34.388006</td>
      <td>26.369286</td>
      <td>3.395453</td>
      <td>5.485167</td>
      <td>6.931863</td>
      <td>0.093706</td>
      <td>19.792547</td>
      <td>38.539763</td>
      <td>25.926938</td>
      <td>27.002840</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 15672 columns</p>
</div>




```python
#### Step3. Generate the inputs and save them
df_inputs = df_cancer_type.join(df_tpm_conceptor)
df_inputs.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cancer_type</th>
      <th>A1BG</th>
      <th>A1CF</th>
      <th>A2M</th>
      <th>A2ML1</th>
      <th>A4GALT</th>
      <th>A4GNT</th>
      <th>AAAS</th>
      <th>AACS</th>
      <th>AADAC</th>
      <th>...</th>
      <th>ZWILCH</th>
      <th>ZWINT</th>
      <th>ZXDA</th>
      <th>ZXDB</th>
      <th>ZXDC</th>
      <th>ZYG11A</th>
      <th>ZYG11B</th>
      <th>ZYX</th>
      <th>ZZEF1</th>
      <th>ZZZ3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ERR2208944</th>
      <td>25</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>859.203620</td>
      <td>73.019466</td>
      <td>11.942279</td>
      <td>1.947147</td>
      <td>86.527503</td>
      <td>9.236956</td>
      <td>3.918524</td>
      <td>...</td>
      <td>16.725820</td>
      <td>14.012390</td>
      <td>7.255890</td>
      <td>6.984732</td>
      <td>16.469669</td>
      <td>0.879873</td>
      <td>21.106355</td>
      <td>89.920944</td>
      <td>47.520979</td>
      <td>21.534480</td>
    </tr>
    <tr>
      <th>ERR2208928</th>
      <td>25</td>
      <td>0.038627</td>
      <td>0.032171</td>
      <td>881.830260</td>
      <td>7.533515</td>
      <td>12.650118</td>
      <td>2.778540</td>
      <td>95.158662</td>
      <td>10.227978</td>
      <td>1.798357</td>
      <td>...</td>
      <td>34.613126</td>
      <td>42.500215</td>
      <td>12.806729</td>
      <td>12.317719</td>
      <td>18.357604</td>
      <td>0.526526</td>
      <td>35.163162</td>
      <td>60.709750</td>
      <td>52.413439</td>
      <td>25.859226</td>
    </tr>
    <tr>
      <th>ERR2208949</th>
      <td>25</td>
      <td>0.030213</td>
      <td>0.012581</td>
      <td>504.984491</td>
      <td>50.836895</td>
      <td>5.900676</td>
      <td>0.611231</td>
      <td>106.174319</td>
      <td>8.090318</td>
      <td>4.960132</td>
      <td>...</td>
      <td>12.677251</td>
      <td>19.670726</td>
      <td>12.836934</td>
      <td>9.511444</td>
      <td>9.528438</td>
      <td>0.154435</td>
      <td>23.424874</td>
      <td>69.710920</td>
      <td>40.638326</td>
      <td>25.926391</td>
    </tr>
    <tr>
      <th>ERR2208900</th>
      <td>25</td>
      <td>0.143579</td>
      <td>0.023916</td>
      <td>1940.416805</td>
      <td>0.182940</td>
      <td>4.014771</td>
      <td>0.813332</td>
      <td>46.225429</td>
      <td>2.235265</td>
      <td>0.042219</td>
      <td>...</td>
      <td>24.998045</td>
      <td>21.527292</td>
      <td>8.889713</td>
      <td>7.719297</td>
      <td>12.813737</td>
      <td>0.670318</td>
      <td>18.452489</td>
      <td>56.563242</td>
      <td>36.542147</td>
      <td>32.326600</td>
    </tr>
    <tr>
      <th>ERR2208922</th>
      <td>25</td>
      <td>0.109992</td>
      <td>0.286277</td>
      <td>1534.682495</td>
      <td>0.860539</td>
      <td>11.887174</td>
      <td>0.061813</td>
      <td>83.537696</td>
      <td>5.875654</td>
      <td>1.482365</td>
      <td>...</td>
      <td>34.388006</td>
      <td>26.369286</td>
      <td>3.395453</td>
      <td>5.485167</td>
      <td>6.931863</td>
      <td>0.093706</td>
      <td>19.792547</td>
      <td>38.539763</td>
      <td>25.926938</td>
      <td>27.002840</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 15673 columns</p>
</div>






```python
df_inputs.to_csv('./data/conceptor_inputs.csv')
```
