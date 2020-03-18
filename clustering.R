install.packages('pheatmap')
install.packages('dendextend')
install.packages('eulerr')

library(tidyverse)
library(pheatmap)
library(dendextend)
library(DESeq2)
library(openxlsx)
library(eulerr)


# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Import in data and subset it for how many genes you want to see.
data <- read.xlsx("DESeq2_exports/caf_age10_LRT_export.xlsx", rowNames = TRUE)

# Recasting character matrices to numeric (if necessary - VERY ANNOYING)
#data <- data.matrix(data, rownames.force = NA)
#data <- data.frame(data)

#Manual Reordering of Columns (if necessary)
reordered_index <- c(
                     # grep("4_wt_APOE3_wt_f", names(data), ignore.case = T),
                     # grep("4_wt_APOE3_wt_m", names(data), ignore.case = T),
                     # grep("4_wt_APOE3_het_f", names(data), ignore.case = T),
                     # grep("4_wt_APOE3_het_m", names(data), ignore.case = T),
                     # grep("4_wt_APOE4_wt_f", names(data), ignore.case = T),
                     # grep("4_wt_APOE4_wt_m", names(data), ignore.case = T),
                     # grep("4_wt_APOE4_het_f", names(data), ignore.case = T),
                     # grep("4_wt_APOE4_het_m", names(data), ignore.case = T),
                     # grep("4_APP_APOE3_wt_f", names(data), ignore.case = T),
                     # grep("4_APP_APOE3_wt_m", names(data), ignore.case = T),
                     # grep("4_APP_APOE3_het_f", names(data), ignore.case = T),
                     # grep("4_APP_APOE3_het_m", names(data), ignore.case = T),
                     # grep("4_APP_APOE4_wt_f", names(data), ignore.case = T),
                     # grep("4_APP_APOE4_wt_m", names(data), ignore.case = T),
                     # grep("4_APP_APOE4_het_f", names(data), ignore.case = T),
                     # grep("4_APP_APOE4_het_m", names(data), ignore.case = T),
                     grep("10_wt_APOE3_wt_f", names(data), ignore.case = T),
                     grep("10_wt_APOE3_wt_m", names(data), ignore.case = T),
                     grep("10_wt_APOE3_het_f", names(data), ignore.case = T),
                     grep("10_wt_APOE3_het_m", names(data), ignore.case = T),
                     grep("10_wt_APOE4_wt_f", names(data), ignore.case = T),
                     grep("10_wt_APOE4_wt_m", names(data), ignore.case = T),
                     grep("10_wt_APOE4_het_f", names(data), ignore.case = T),
                     grep("10_wt_APOE4_het_m", names(data), ignore.case = T),
                     grep("10_APP_APOE3_wt_f", names(data), ignore.case = T),
                     grep("10_APP_APOE3_wt_m", names(data), ignore.case = T),
                     grep("10_APP_APOE3_het_f", names(data), ignore.case = T),
                     grep("10_APP_APOE3_het_m", names(data), ignore.case = T),
                     grep("10_APP_APOE4_wt_f", names(data), ignore.case = T),
                     grep("10_APP_APOE4_wt_m", names(data), ignore.case = T),
                     grep("10_APP_APOE4_het_f", names(data), ignore.case = T),
                     grep("10_APP_APOE4_het_m", names(data), ignore.case = T)
)

reordered_data <- data %>%
  select(reordered_index)

# Data subsetting and gene filtering (if necessary)
unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps'), collapse = '|')

data_final <- reordered_data %>%
  rownames_to_column('gene') %>%
  filter(!str_detect(gene, unwanted_genes)) %>%
  column_to_rownames('gene')
data_subset <- data_final %>% head(n = 100)

# Data grouping by means (if desired)
data_grouped <- data.frame(
  # wt3wt_4 = rowMeans(select(data_final, starts_with(("4_wt_APOE3_wt")))),
  # wt3het_4 = rowMeans(select(data_final, starts_with(("4_wt_APOE3_het")))),
  # wt4wt_4 = rowMeans(select(data_final, starts_with(("4_wt_APOE4_wt")))),
  # wt4het_4 = rowMeans(select(data_final, starts_with(("4_wt_APOE4_het")))),
  # APP3wt_4 = rowMeans(select(data_final, starts_with(("4_APP_APOE3_wt")))),
  # APP3het_4 = rowMeans(select(data_final, starts_with(("4_APP_APOE3_het")))),
  # APP4wt_4 = rowMeans(select(data_final, starts_with(("4_APP_APOE4_wt")))),
  # APP4het_4 = rowMeans(select(data_final, starts_with(("4_APP_APOE4_het")))),
  wt3wt_10 = rowMeans(select(data_final, starts_with(("10_wt_APOE3_wt")))),
  wt3het_10 = rowMeans(select(data_final, starts_with(("10_wt_APOE3_het")))),
  wt4wt_10 = rowMeans(select(data_final, starts_with(("10_wt_APOE4_wt")))),
  wt4het_10 = rowMeans(select(data_final, starts_with(("10_wt_APOE4_het")))),
  APP3wt_10 = rowMeans(select(data_final, starts_with(("10_APP_APOE3_wt")))),
  APP3het_10 = rowMeans(select(data_final, starts_with(("10_APP_APOE3_het")))),
  APP4wt_10 = rowMeans(select(data_final, starts_with(("10_APP_APOE4_wt")))),
  APP4het_10 = rowMeans(select(data_final, starts_with(("10_APP_APOE4_het"))))
)
data_grouped_subset <- data_grouped %>% head(n=100)

# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_subset_norm <- t(apply(data_subset, 1, cal_z_score))


## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Organize your dendrograms and cluster groupings
#my_hclust_gene <- hclust(dist(data_subset), method = 'complete')
#my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
#my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))

#my_sample_col <- data.frame(sample = rep(c('male', 'female'), c))

# Plot gene dendrogram to see what it looks like
#as.dendrogram(my_hclust_gene) %>%
#  plot(horiz = TRUE)


# Create heatmap
phm_full <- pheatmap(data_norm,
                color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                breaks = seq(-2, 2, by = 0.1),
                kmeans_k = NA,
                cluster_rows = T,
                cutree_row = 3,
                cluster_cols = F,
                #cutree_cols = 3,
                #gaps_col = c(22,45,63,84,104,127,144),  # raw
                #gaps_col = c(20,43,55,76,94,115,130),  # raw filtered
                #gaps_col = c(11,22,31,39,50,62,72),    # O
                #gaps_col = c(11,23,32,45,54,65,72),    # P
                #gaps_col = c(8,17,28,36,45,51,59),     # APP
                #gaps_col = c(12,20,32,43,53,60,70),    # age10
                gaps_col = c(11,18,30,39,48,54,61),    # age10 filtered
                legend = TRUE,
                show_rownames = F,
                treeheight_col = 0,
                treeheight_row = 0,
                scale = 'row'
)

phm <- pheatmap(data_subset_norm,
                color = colorRampPalette(c("navy", "white", "firebrick3"))(61),
                breaks = seq(-3, 3, by = 0.1),
                kmeans_k = NA,
                cluster_rows = T,
                cutree_row = 1,
                cluster_cols = F,
                cutree_cols = 3,
                #gaps_col = c(22,45,63,84,104,127,144),  # raw
                #gaps_col = c(20,43,55,76,94,115,130),  # raw filtered
                #gaps_col = c(11,22,31,39,50,62,72),    # O
                #gaps_col = c(11,23,32,45,54,65,72),    # P
                #gaps_col = c(8,17,28,36,45,51,59),     # APP
                #gaps_col = c(12,20,32,43,53,60,70),    # age10
                gaps_col = c(11,18,30,39,48,54,61),    # age10 filtered
                legend = TRUE,
                show_rownames = T,
                treeheight_col = 0,
                treeheight_row = 50,
                scale = 'row'
)
                
# Extract cluster genes (if desired)
# Build cluster tables from dendrogram based on k
gene_clusters <- data.frame(sort(cutree(phm$tree_row, k=2))) %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

gene_clusters_br <- gene_clusters %>% rename(cluster = sort.cutree.phm.tree_row..k...2..)

# Cluster Comparison
gene_clusters_scXbr <- as_tibble(full_join(gene_clusters_br, gene_clusters_sc, by = 'gene'))
gene_clusters_scXbr <- gene_clusters_scXbr %>% 
  mutate(cluster.x = recode(cluster.x, '1' = 'down', '2' = 'up')) %>%
  mutate(cluster.y = recode(cluster.y, '1' = 'down', '2' = 'up'))
shared_DEgenes_scXbr <- gene_clusters_scXbr %>%
  filter(!is.na(cluster.x) & !is.na(cluster.y))
  
vd <- euler(c('Brain' = 142, 'Spinal Cord' = 160, 'Shared' = 9))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Build Cluster plots
cluster_rlog <- rld_mat[gene_clusters$gene, ]

# Build cluster plot by variable (time parameter)
cluster_group <- degPatterns(cluster_rlog, metadata = sod_noPBS_data, time = 'treatment', col = 'tissue_type', reduce = T)

# Write to file
write_tsv(gene_clusters, 'results/female_condition_DE_clusters.tsv')
write_tsv(pain_upreg_genes_male, 'results/pain_upreg_genes_male')
write_tsv(gene_clusters_scXbr, 'results/xCondition_noPBS_gene_clusters_scXbr.tsv')
write_tsv(shared_DEgenes_scXbr, 'results/xCondition_noPBS_shared_DEGenes_scXbr.tsv')


