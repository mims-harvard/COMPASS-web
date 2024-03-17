theme: jekyll-theme-minimal

## Welcome to conceptor pages

## **conceptor 101**
#### **The documents and instructions for the immnuethrepy response prediction model, conceptor**


**Date: 2024-03-17**



# **Content**
[Input data preparation](https://zitniklab.hms.harvard.edu/conceptor-101/00_prepare_input_data.html)

[1. Make predictions](#make_prediction)



# **Input data preparation**

## Conceptor Input Requirement

Please note that the input for Conceptor should be mRNA's TPM (**Transcripts Per Million**) expression values, not raw counts or other forms of mRNA expression values. TPM calculation is similar to FPKM but differs in the normalization process. In TPM, all transcripts are normalized for length first. Then, instead of using the total overall read count for size normalization, the sum of the length-normalized transcript values is used as a size indicator.

Please fell free to contact me if you have any questions on this.

## Data Processing Recommendation

If your data is in raw sequence format (FASTQ) or as raw counts, we recommend processing it using the following bioinformatics pipeline. This recommendation is based on the fact that our pretrained TCGA (The Cancer Genome Atlas) data was processed using this pipeline, and using the same pipeline for your input data may yield better results.


Please fell free to contact me if you have any questions on this.


### 1. [mRNA-Seq Alignment Workflow](https://github.com/mims-harvard/conceptor-101/tree/main/docs/mRNA_pipeline)

The RNA-Seq Alignment Workflow follows these steps:

`fastqc/multiqc  --> fastp ---> STAR2 align (two-pass method)`

For more information, please refer to [GDC mRNA expression pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/).


![pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/images/RNA-Seq-DR32_Image.png)

#### Specific Process

1. Begin with quality control using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiqc](http://multiqc.info).
2. Proceed with data preprocessing using [fastp](https://github.com/OpenGene/fastp) to clean raw sequencing data and improve quality.
3. Finally, align the RNA-seq reads to a reference genome using STAR version 2.7.5c, which maps RNA-seq reads to the [reference genome](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files). While custom index files can be created, we use the reference genome files downloaded from GDC. The link for the specific reference genome file (star-2.7.5c_GRCh38.d1.vd1_gencode.v36.tgz) is available [here](https://api.gdc.cancer.gov/data/c0008693-0583-4eac-bd5c-583070763893).



### 2. [Converting mRNA Raw Counts to TPM](./mRNA_pipeline)

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

# convert count to TPM based the gtf file
gtf_path = "./data/gencode.v36.annotation.gtf.gz"
tpm = TPM(gtf_path).set_output(transform="pandas")

# example of the raw counts
df_counts = pd.read_csv('./data/toy_raw_counts.csv', index_col=0)
df_counts.head()

# example of the TPM values
df_tpm = tpm.fit_transform(df_counts)
df_tpm.to_csv('./data/toy_tpm.csv')
df_tpm.head()
```


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






# **make_prediction**

## **2.1 Item1**
## **2.2 Item2**
