library("SummarizedExperiment")
library("dplyr")

acLDL_se = readRDS("results/SummarizedExperiments/acLDL_salmon_Ensembl_87.rds")
other_se = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
expressed_genes = unique(c(rowData(acLDL_se)$gene_id, rowData(other_se)$gene_id))

#Import QTL counts
ctrl_vs_acldl = readr::read_delim("results/DE/Ctrl_vs_AcLDL_DESeq2_fold_change.txt", 
                                  delim = "\t", col_types = "ccdddddd") %>%
  dplyr::filter(gene_id %in% expressed_genes)
naive_vs_IFNg = readr::read_delim("results/DE/naive_vs_IFNg_DESeq2_fold_change.txt.gz", 
                                  delim = "\t", col_types = "ccdddddd") %>%
  dplyr::filter(gene_id %in% expressed_genes)
naive_vs_SL = readr::read_delim("results/DE/naive_vs_Salmonella_DESeq2_fold_change.txt.gz", 
                                  delim = "\t", col_types = "ccdddddd") %>%
  dplyr::filter(gene_id %in% expressed_genes)
naive_vs_IFNg_SL = readr::read_delim("results/DE/naive_vs_IFNg+Salmonella_DESeq2_fold_change.txt.gz", 
                                  delim = "\t", col_types = "ccdddddd") %>%
  dplyr::filter(gene_id %in% expressed_genes)

de_counts = list(AcLDL = nrow(dplyr::filter(ctrl_vs_acldl, padj < 0.01, abs(log2FoldChange) > 1)),
     SL1344 = nrow(dplyr::filter(naive_vs_SL, padj < 0.01, abs(log2FoldChange) > 1)),
     IFNg = nrow(dplyr::filter(naive_vs_IFNg, padj < 0.01, abs(log2FoldChange) > 1)),
     IFNg_SL1344 = nrow(dplyr::filter(naive_vs_IFNg_SL, padj < 0.01, abs(log2FoldChange) > 1)))
saveRDS(de_counts, "results/DE/de_counts.rds")

