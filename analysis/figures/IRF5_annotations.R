library("wiggleplotr")
library("GenomicRanges")
library("dplyr")
library("devtools")
load_all("../txrevise/")
library("ggplot2")

#Import IRF5 transcripts
IRF5_data = readRDS("../txrevise/data/IRF5.rds")
plotting_annotations = dplyr::select(IRF5_data$metadata, ensembl_transcript_id, ensembl_gene_id, external_gene_name, strand) %>% 
  dplyr::rename(transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)

irf5_tx = wiggleplotr::plotTranscripts(IRF5_data$exons, IRF5_data$cdss, rescale_introns = TRUE, transcript_label = FALSE)
ggsave("results/figures/IRF5_annotated_transcripts.pdf", plot = irf5_tx, width = 3, height = 4)

irf5_tx = wiggleplotr::plotTranscripts(IRF5_data$exons, IRF5_data$cdss, plotting_annotations, rescale_introns = TRUE)
ggsave("results/figures/IRF5_annotated_transcripts_labeled.pdf", plot = irf5_tx, width = 5, height = 6)

#Extend truncated transcripts
gene_extended_tx = txrevise::extendTranscriptsPerGene(IRF5_data$metadata, IRF5_data$exons, IRF5_data$cdss)
gene_data_ext = txrevise::replaceExtendedTranscripts(IRF5_data, gene_extended_tx)
irf5_extended_tx = wiggleplotr::plotTranscripts(gene_data_ext$exons, gene_data_ext$cdss, plotting_annotations, 
                             rescale_introns = TRUE)
ggsave("results/figures/IRF5_extended_transcripts.pdf", plot = irf5_extended_tx, width = 5, height = 6)

