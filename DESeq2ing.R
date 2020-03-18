# References: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

## Gene-level differential expression analysis using DESeq2
## !requireNamespace is used to check whether the package is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

##biocManager is a package, and BiocManager::install is used to install packages
BiocManager::install("BiocGenerics")
BiocManager::install("DESeq2")
BiocManager::install("tidyverse")
BiocManager::install("DEGreport")
## tximport is to import transcript-level abundances and estimated counts for gene-level analysis packages
BiocManager::install("tximport")
BiocManager::install("ggplot2")
BiocManager::install("biomaRt")
BiocManager::install("apeglm")
BiocManager::install("limma")

library(BiocGenerics)
library(tximport)
library(S4Vectors)
library(DESeq2)
library(openxlsx)
library(biomaRt)
library(tidyverse)
library(DEGreport)
library(ggplot2)
library(apeglm)
library(limma)


##---------------------------------FILE PREPARATION----------------------------------##

## useMart is a function used to connect to the selected BioMart database and dataset
ensMart <- useMart("ensembl")
ensembl_ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
## getBM is a function used to retrieve information from the BioMart database
ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), mart = ensembl_ms_mart) %>%
  as_tibble()

# Build tx2gene table from ensembl_df
tx2gene <- ensembl_df %>%
  select(ensembl_transcript_id_version, external_gene_name) %>%
  rename(Transcript = ensembl_transcript_id_version, Gene = external_gene_name) %>%
  data.frame()

## Create a metadata map using your metadata file.
metadata <- read.xlsx("meta_app.xlsx", sheet = 1)
sample_data <- data.frame(position = metadata$Position, id = metadata$Sample_ID, gender = metadata$Gender, age = metadata$Sacrifice, genotype = metadata$APP, APOE = metadata$APOE, CX3 = metadata$CX3CR1creERT2, plate = metadata$Plate)
sample_data$age <- factor(sample_data$age)

# Outlier filtering  and subsetting(if necessary)
#outliers = c('')
outliers = c('Fp10', 'A08', 'Bp03', 'Bp09', 'C01', 'D04', 'F11', 'Dp10', 'B02', 'C08', 'C11', 'Ep02', 'A01', 'Ap05', 'E01') 
sample_data <- sample_data %>%            # Filter sample_data sheet
  filter(!position %in% outliers)

age10_data <- sample_data %>% filter(age == '10')
app_data <- sample_data %>% filter(genotype == 'APP')
o_data <- sample_data %>% filter(plate == 'O')
p_data <- sample_data %>% filter(plate == 'P')

## List all directories containing data
#*** All the gene quantification files are in the folder "quant_files", make sure metafile sample order matches order of quant files.
all_files <- list.files("./quant_files", full.names = T, pattern=NULL, all.files=FALSE)
quant_files <- file.path(all_files, "quant.sf")

position_list <- paste0(sample_data$position, collapse = '|')
sample_files <- grep(position_list, quant_files, value = TRUE)


## QC check - Order of metafile samples must match quant file order!  Fix metadata file if not!
all(file.exists(quant_files))
sample_files
sample_data

# Import Transcript quantification files
txi<-tximport(files=sample_files,type="salmon", tx2gene = tx2gene)

## Look at the counts and set the counts as a matrix file.
txicounts<-as.matrix(txi$counts)

## Write the counts to file, round them up, and convert to data.frame and remove txicounts file.
count_data <- txi$counts %>%
  round() %>%
  data.frame()
rm(txicounts)


#-----------------------------------ANALYSIS BEGINS HERE---------------------------------#


## Create DESeq2Dataset object with variables of interest.
dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design = ~ plate + gender + age + genotype + APOE + CX3 + gender:genotype + gender:APOE + genotype:APOE + genotype:CX3) 

#prefiltering on minimum of 10 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# Set up groups for the analysis.
#dds$Condition <- relevel(dds$Condition, ref = "Control")
#dds$Region <- relevel(dds$Region, ref = "Ipsi/lum/sc")
#dds$Region <- factor(dds$Region, levels = c("PAG/RVM", "PTN/SI", "PFC/HP/AMY", "Ipsi/lum/sc"))
#dds$Gender <- factor(dds$Gender, levels = c("m", "f"))

dds$group <- factor(paste0(dds$tissue_type, "_", dds$treatment))
design(dds) <- ~group
dds$group

## Run DESeq analysis to gather differential expression results
# Run DESeq (Wald)
dds_Wald <- DESeq(dds)
# Or run with DESeq (LRT)
dds_LRT <- DESeq(dds, test = 'LRT', reduced = ~ plate + gender + age)

# View names of estimated effects
resultsNames(dds_Wald)
resultsNames(dds_LRT)

# Create table of effects showing log2fold, standard error, test stats, and p-vals
res_Wald <- results(dds_Wald)
res_LRT <- results(dds_LRT, pAdjustMethod = 'BH')
#res_sod_PBS <- results(dds, name = 'treatment_SODMM_vs_PBS')

summary(res_Wald)
summary(res_LRT)

## View Total number of normalized counts per sample
dds_LRT <- estimateSizeFactors(dds_LRT)
DS_norm_counts <- counts(dds_LRT, normalized = TRUE)

#colSums(DS_norm_counts)

## Plot dispersion estimates
plotDispEsts(dds_LRT)

# Transform counts for data visualization
rld <- vst(dds_LRT, blind=TRUE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch, design=~Condition)

# Plot PCA 
plotPCA(rld, intgroup=c('genotype', 'plate', 'CX3'))

# Add nametags
z <- plotPCA(rld, intgroup=c('genotype', 'CX3'))
z + geom_label(aes(label = age10_data$position))

# Build significant gene table and extract list of sorted DE genes.
sig_res <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%    # makes a 'gene' column using column1
  as_tibble %>%
  filter(padj < 0.1) %>%
  arrange(padj)
sorted_DEGenes <- sig_res$gene

# Export significant gene count data
name_list <- c('gene', (paste0(dds_LRT$age, "_", dds_LRT$genotype, "_", dds_LRT$APOE, "_", dds_LRT$CX3, "_", dds_LRT$gender, "_", dds_LRT$position)))

sig_app_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  slice(match(sorted_DEGenes, gene))
  

# Write DE genes to file
write.xlsx(sig_app_export, 'DESeq2_exports/caf_age10_LRT_export.xlsx')

## Define contrasts, extract results table, and shrink the log2 fold changes

#contrast_1 <- c("Group", "Control_PAG.RVM_f", "Control_Ipsi.lum.sc_f")
#resPFC_HP_AMYm_vs_Ipsi_lum_scf <- results(dds, contrast=contrast_1, alpha = 0.05)
#resPFC_HP_AMYm_vs_Ipsi_lum_scf_ordered <- resPFC_HP_AMYm_vs_Ipsi_lum_scf[order(resPFC_HP_AMYm_vs_Ipsi_lum_scf$padj),]


#summaries
#summary(resPFC_HP_AMYm_vs_Ipsi_lum_scf_ordered)
#plotMA(resPFC_HP_AMYm_vs_Ipsi_lum_scf_ordered)

#res1 <- lfcShrink(dds, coef="group_PFC_HP_AMYm_vs_Ipsi_lum_scf", type = "apeglm")
#plotMA(res1)

#plotcounts
#setwd("/Users/butovskylab/Desktop/pain\ project/results")
#plotCounts(dds, gene=which.min(resPFC_HP_AMYm_vs_Ipsi_lum_scf$padj), intgroup="Genotype")
#mcols(resPFC_HP_AMYm_vs_Ipsi_lum_scf)$description
#write.csv(as.data.frame(resPFC_HP_AMYm_vs_Ipsi_lum_scf_ordered), 
#          file="condition_PFC_HP_AMYm_vs_Ipsi_lum_scf_results0.05.csv")
#save(resPFC_HP_AMYm_vs_Ipsi_lum_scf, dds, sampleTable, file="DEPFC_HP_AMYm_vs_Ipsi_lum_scf.RData")


# Troubleshooting
# Column Names
check_cols <- c('gene', colnames(data_final))

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  select(c(1, reordered_index+1)) %>% 
  as_tibble %>%
  `colnames<-`(check_cols)

# Results
check_res <- res_LRT %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  filter(!str_detect(gene, unwanted_genes))

# TPM File
my_genes <- rownames(txi$abundance)
my_genes_ann <- ensembl_df[match(my_genes,ensembl_df$external_gene_name),]
app_TPM <- cbind(my_genes_ann$external_gene_name,txi$abundance) %>%
  as.data.frame() %>%
  select(c(1, reordered_index+1)) %>%
  as_tibble() %>%
  `colnames<-`(check_cols)

write.xlsx(check_res, file="results/caf_age10_LRT_results_filtered.xlsx")
write.xlsx(check_counts, file='results/caf_age10_LRT_DS_counts.xlsx')
write.xlsx(app_TPM, file="results/caf_age10_LRT_TPM.xlsx")



