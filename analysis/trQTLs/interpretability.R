library("dplyr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")

#Count alternative transcripts
revised_count = rowData(se_revised) %>% tbl_df2() %>% dplyr::group_by(group_id) %>% 
  dplyr::summarise(transcript_count = length(transcript_id)) %>% 
  dplyr::mutate(quant = "reviseAnnotations")
ensembl_count = rowData(se_ensembl) %>% tbl_df2() %>% dplyr::group_by(gene_id) %>% 
  dplyr::summarise(transcript_count = length(transcript_id)) %>% 
  dplyr::mutate(quant = "Ensembl_87") %>%
  dplyr::filter(transcript_count > 1)
leafcutter_count = rowData(se_leafcutter) %>% tbl_df2() %>% dplyr::group_by(gene_id) %>% 
  dplyr::summarise(transcript_count = length(transcript_id)) %>% 
  dplyr::mutate(quant = "leafcutter")

#Merge counts
tx_count = bind_rows(revised_count, ensembl_count, leafcutter_count) %>%
  dplyr::select(transcript_count, quant) %>% 
  dplyr::left_join(phenotypeFriendlyNames())
medians = dplyr::group_by(tx_count, phenotype) %>% dplyr::summarise(median_tx = median(transcript_count))

#Make histogram
histogram = ggplot(tx_count, aes(x = transcript_count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim = c(0,30)) + 
  facet_grid(phenotype~., scales = "free_y") +
  theme_light() + 
  theme(legend.position = "none") +
  xlab("Number of alternative features")
ggsave("results/figures/trQTL_interpretability.pdf", histogram, width = 4, height = 4)
