library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("data.table")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr")

#Import SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")

#Import GRangesList
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")

#Extract transcript metadata
tx_metadata = rowData(se_revised) %>% tbl_df2() %>% 
  dplyr::filter(gene_name == "IRF5") %>% 
  tidyr::separate(group_id, c("gene","group","position"), "\\.") %>%
  dplyr::select(-gene) %>%
  dplyr::filter(group == "grp_1")

#Upstream
txs = dplyr::filter(tx_metadata, position == "upstream")
tx_ranges = revised_granges[txs$transcript_id]
plotTranscripts(tx_ranges)

#Contained
txs = dplyr::filter(tx_metadata, position == "contained")
tx_ranges = revised_granges[txs$transcript_id]
plotTranscripts(tx_ranges)

#Downstream
txs = dplyr::filter(tx_metadata, position == "downstream")
tx_ranges = revised_granges[txs$transcript_id]
plotTranscripts(tx_ranges)


