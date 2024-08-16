#Created on Wed Aug 16 10:38:40 2023
#@author: Wanxiang Shen

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data") ## adapte to new version GDC
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("maftools")
# BiocManager::install("mclust")
# BiocManager::install("BiocOncoTK")
# BiocManager::install("GSVA")



library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)
library(BiocOncoTK)

dirpath <- '/n/data1/hms/dbmi/zitnik/lab/users/was966/TCGA/'
directory = paste0(dirpath, 'GDC_v37')


if (!dir.exists(dirpath)){
  dir.create(dirpath)
}

if (!dir.exists(directory)){
  dir.create(directory)
}

setwd(dirpath)

#version_https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/#data-release-370
tcga_list = c('TCGA-BRCA','TCGA-GBM','TCGA-OV','TCGA-LUAD','TCGA-UCEC','TCGA-KIRC',
              'TCGA-HNSC','TCGA-LGG','TCGA-THCA','TCGA-LUSC','TCGA-PRAD','TCGA-SKCM',
              'TCGA-COAD','TCGA-STAD','TCGA-BLCA','TCGA-LIHC','TCGA-CESC','TCGA-KIRP',
              'TCGA-SARC','TCGA-LAML','TCGA-ESCA','TCGA-PAAD','TCGA-PCPG','TCGA-READ',
              'TCGA-TGCT','TCGA-THYM','TCGA-KICH','TCGA-ACC','TCGA-MESO','TCGA-UVM',
              'TCGA-DLBC','TCGA-UCS','TCGA-CHOL')


#save SMI table
data(MSIsensor.10k)
data(patient_to_tumor_code)
names(MSIsensor.10k)[1] = "patient_barcode"
write.table(MSIsensor.10k, file = paste0(directory, '/SMI_table.txt'), row.names = F, sep = "\t", quote = F)



for(i in tcga_list){

    ## expression ###
    query <- GDCquery(project = i,
                    data.category = 'Transcriptome Profiling',
                    data.type = 'Gene Expression Quantification',
                    workflow.type = 'STAR - Counts',) #, "Primary Tumor" #sample.type = c("Metastatic")

    if(is.null(query) != TRUE){
        GDCdownload(query,directory = directory)
        data <- GDCprepare(query, directory = directory)
        #use assayNames(data) to select the tpm_unstrand column

        data_counts <- assay(data,i = 'unstranded')#
        data_tpm <- assay(data,i = 'tpm_unstrand')#
    
        rowdata <- rowData(data) ## gene table
        coldata <- colData(data) ## patient table

        #treatments <- coldata$treatments
        #primary_site <- coldata$primary_site
        #primary_site_df = as.data.frame(do.call(rbind, primary_site))
        #treatments_df = as.data.frame(do.call(rbind, treatments))
    
        # remove \n in treatments col of coldata
        coldata$treatments <- gsub("\\n","",coldata$treatments )
        coldata$primary_site <- gsub("\\n","",coldata$primary_site )
    
        write.csv(data_counts, file = paste0(directory, '/', i,'/rnaSeq_counts_matrix.csv'))
        write.csv(data_tpm, file = paste0(directory, '/', i, '/rnaSeq_tpm_matrix.csv'))
        write.table(rowdata, file = paste0(directory, '/', i, '/rnaSeq_gene_table.txt'), row.names = F, sep = "\t", quote = F)
        write.table(coldata, file = paste0(directory, '/', i, '/rnaSeq_sample_table.txt'), row.names = F, sep = "\t", quote = F)
        
    }

    ### simple nucleotide variation ###
    try_res <- try(querySNV <- GDCquery(project = i, 
                     data.category = "Simple Nucleotide Variation", # Simple nucleotide variation if legacy
                     data.type = "Masked Somatic Mutation",
                     access = "open", 
                     data.format = "MAF", #sample.type = "Primary Tumor"
                     experimental.strategy = c("WXS")),  silent=TRUE)

    if (! inherits(try_res, 'try-error')) { 
        if(is.null(querySNV) != TRUE){
            GDCdownload(querySNV, directory = directory)
            SNVdata <- GDCprepare(querySNV, directory = directory)
            #write.csv(SNVdata, file = paste0(directory, '/', i, '/SNV_table.csv')
            SNVdata$Tumor_Sample_Barcode = substr(SNVdata$Tumor_Sample_Barcode, 1, 12)
            SNVdata$t_vaf = SNVdata$t_alt_count/SNVdata$t_depth * 100
            write.table(SNVdata, file = paste0(directory, '/', i, '/SNV_table.txt'), row.names = F, sep = "\t", quote = F)
        }
        }

    
    ### Copy Number Variation ###
    try_res <- try(queryCNV <- GDCquery(project = i, 
                         data.category = "Copy Number Variation",
                         data.type = "Masked Copy Number Segment", 
                         access = "open"),  silent=TRUE)
    if (! inherits(try_res, 'try-error')) {
        if(is.null(queryCNV) != TRUE){
            GDCdownload(queryCNV, directory = directory)
            CNVdata <- GDCprepare(queryCNV, directory = directory)
            CNVdata$Sample = substr(CNVdata$Sample, 1, 12)
            seg <- CNVdata[, c(7, 2:6)]
            seg = data.table::data.table(seg)
            head(seg)
            write.table(seg, file = paste0(directory, '/', i, '/SCNV_table.txt'), row.names = F, sep = "\t", quote = F)
            }}}


    
    ### clinical ###
    try_res <- try(clinical_query <- GDCquery(project = i,
                    data.category = 'Clinical',
                    data.format = 'bcr xml'),  silent=TRUE)
    if (! inherits(try_res, 'try-error')) {
        if(is.null(clinical_query) != TRUE){    
            GDCdownload(clinical_query, directory = directory)
            ## patient information ##
            patient <- GDCprepare_clinic(clinical_query, clinical.info = "patient", directory = directory)
            write.table(patient, file = paste0(directory, '/', i, '/clinical_patient_table.txt'), row.names = F, sep = "\t", quote = F)
            ## drug information ##
            drug <- GDCprepare_clinic(clinical_query, clinical.info = "drug", directory = directory)
            write.table(drug, file = paste0(directory, '/', i, '/clinical_drug_table.txt'), row.names = F, sep = "\t", quote = F)
        }}


    ## Proteome Profiling ##
    try_res <- try(prot_query <- GDCquery(
                                project = i,
                                data.category = "Proteome Profiling",
                                data.type = "Protein Expression Quantification"),  silent=TRUE)
    
    if (! inherits(try_res, 'try-error')) {
        if(is.null(prot_query) != TRUE){
            GDCdownload(prot_query, directory = directory) 
            #### remove duplicates by cases
            prot_query.1=prot_query
            tmp=prot_query.1$results[[1]]
            tmp=tmp[which(!duplicated(tmp$cases)),]
            prot_query.1$results[[1]]=tmp
            proteome <- GDCprepare(prot_query.1, directory = directory)
            write.table(proteome, file = paste0(directory, '/', i, '/proteome_table.txt'), row.names = F, sep = "\t", quote = F)}}
}

