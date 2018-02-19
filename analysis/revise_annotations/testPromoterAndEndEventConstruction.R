library("devtools")
library("dplyr")
library("data.table")
library("SummarizedExperiment")
library("ggplot2")
load_all("../txrevise/")
load_all("../seqUtils/")
library("wiggleplotr")

#Import revised transcript annotations
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")

#Import gene metadata
revised_gene_metadata = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2()
group_map = dplyr::select(revised_gene_metadata, transcript_id, group_id)


#Make plots for IRF2 transcripts before and after filling in alt exons
#IRF2
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "IRF2") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
original_promoters = plotTranscripts(transcripts)
ggsave("results/figures/IRF2_original_promoters.pdf", plot = original_promoters, width = 5, height = 4)

new_transcripts = fillMissingInternalExons(transcripts, type = "start")
new_promoters = plotTranscripts(new_transcripts)
ggsave("results/figures/IRF2_new_promoters.pdf", plot = new_promoters, width = 5, height = 3)


##### Tests cases for new promoter code ####

#IRF5
#Identify all promoter events
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "IRF5") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "start")
plotTranscripts(new_transcripts)

#KIAA1671
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "KIAA1671") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "start")
plotTranscripts(new_transcripts)

#IRF2
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "IRF2") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "start")
plotTranscripts(new_transcripts)


#CD40
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "CD40") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "start")
plotTranscripts(new_transcripts)

#SLC9B2
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "SLC9B2") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "start")
plotTranscripts(new_transcripts)



#Perform the same analysis for 3'ends

#Identify all promoter events
end_events = dplyr::filter(revised_gene_metadata, gene_name == "IRF5") %>% 
  dplyr::filter(group_id %like% "downstream") %>% 
  dplyr::filter(group_id %like% "grp_2")

transcripts = revised_granges[end_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "end")
plotTranscripts(new_transcripts)


end_events = dplyr::filter(revised_gene_metadata, gene_name == "XRN1") %>% 
  dplyr::filter(group_id %like% "downstream") %>% 
  dplyr::filter(group_id %like% "grp_2")

transcripts = revised_granges[end_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts, type = "end")
plotTranscripts(new_transcripts)

