library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
library("GenomicFeatures")
load_all("analysis/housekeeping/")
load_all("../seqUtils/")
load_all("../txrevise/")

#Import minimal QTL p-values
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
acldl_qtls = readRDS("results/trQTLs/acLDL_trQTL_min_pvalues.rds")

#Import transcript expression
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
tpm_matrix = assays(se_ensembl)$tpms
sample_metadata = as.data.frame(colData(se_ensembl))
gene_metadata = tbl_df2(rowData(se_ensembl))

#Import transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names=TRUE)
cdss = cdsBy(txdb, by = "tx", use.names=TRUE)

#Calculate mean expression in each condition
mean_tpm_per_condition = calculateMean(tpm_matrix, sample_metadata, "condition_name")
mean_tpm_df = dplyr::mutate(mean_tpm_per_condition, transcript_id = rownames(mean_tpm_per_condition)) %>%
  dplyr::left_join(dplyr::select(gene_metadata, gene_id, transcript_id), by = "transcript_id") %>%
  dplyr::transmute(expressed_transcript_id = transcript_id, gene_id, tpm_naive = naive) %>%
  tbl_df2()

#Identify QTLs
naive_qtls = dplyr::filter(salmonella_qtls$Ensembl_87$naive, p_fdr < 0.1)
qtl_transcripts = dplyr::transmute(naive_qtls, gene_id = group_id, qtl_transcript_id = phenotype_id)

#Find the second most highly expressed transcript
expressed_transcripts = dplyr::filter(mean_tpm_df, gene_id %in% qtl_transcripts$gene_id) %>% 
  dplyr::filter(!(expressed_transcript_id %in% qtl_transcripts$qtl_transcript_id)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(gene_id, -tpm_naive) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
all_transcript_pairs = dplyr::left_join(expressed_transcripts, qtl_transcripts)


findAllDiffs <- function(tx1, tx2, exons){
  print(paste(tx1, tx2))
  diff = txrevise::indentifyAddedRemovedRegions(tx1, tx2, exons)
  if(!is.null(diff)){
    all_diffs = dplyr::bind_rows(as.data.frame(elementMetadata(diff[[1]])), as.data.frame(elementMetadata(diff[[2]]))) %>%
      dplyr::summarise(upstream = sum(upstream), downstream = sum(downstream), contained = sum(contained)) %>%
      dplyr::mutate(tx1_id = tx1, tx2_id = tx2) %>%
      dplyr::select(tx1_id, tx2_id, everything())
    return(all_diffs)
  } else{
    return(NULL)
  }
}

#Find all differences between the two transcripts
tx1_list = as.list(all_transcript_pairs$expressed_transcript_id)
tx2_list = as.list(all_transcript_pairs$qtl_transcript_id)
all_differences = purrr::map2(tx1_list, tx2_list, ~findAllDiffs(.x, .y, exons)) %>% purrr::map_df(identity)

#Count the total number of differences
all_diff_df = dplyr::mutate(all_differences, diff_count = rowSums(sign(all_differences[,c(3:5)])))

#Estimate fraction
table(all_diff_df$diff_count)/sum(table(all_diff_df$diff_count))
