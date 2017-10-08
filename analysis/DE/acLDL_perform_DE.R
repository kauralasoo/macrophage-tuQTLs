library("devtools")
library("dplyr")
library("ggplot2")
library("DESeq2")
load_all("../seqUtils/")

#### Import data ####
#Load the raw eQTL dataset
se_featureCounts = readRDS("results/SummarizedExperiments/acLDL_featureCounts.rds")
gene_name_map = dplyr::select(rowData(se_featureCounts) %>% tbl_df2(), gene_id, gene_name)

#### DESeq2 ####
#Use DESeq to identify genes that vary between conditions
design = colData(se_featureCounts) %>% tbl_df2() %>% as.data.frame()
rownames(design) = design$sample_id

dds = DESeq2::DESeqDataSetFromMatrix(assays(se_featureCounts)$counts, design, ~condition_name) 
dds = DESeq2::DESeq(dds, test = "LRT", reduced = ~ 1)
saveRDS(dds, "results/DE/acLDL_DESeq2_condition_name_LRT_results.rds")

#Extract differentially expressed genes in each condition
acLDL_genes = DESeq2::results(dds, contrast=c("condition_name","AcLDL","Ctrl")) %>% 
  tidyDESeq(gene_name_map) 

#Save DE tables to disk
write.table(acLDL_genes, "results/DE/Ctrl_vs_AcLDL_DESeq2_fold_change.txt", sep = "\t", quote = FALSE, row.names = FALSE)

