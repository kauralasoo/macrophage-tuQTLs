library("devtools")
library("dplyr")
library("data.table")
library("SummarizedExperiment")
load_all("../txrevise/")
load_all("../seqUtils/")
library("wiggleplotr")

#Helper functions
calculatePairDiff <- function(tx1_id, tx2_id, exons_list){
  diff = txrevise::indentifyAddedRemovedRegions(tx1_id, tx2_id, exons_list)
  diff_df = elementMetadata(c(diff[[1]], diff[[2]])) %>% as.data.frame() %>% 
    colSums() %>% as.data.frame() %>% t() %>% dplyr::tbl_df() %>%
    dplyr::mutate(tx1 = tx1_id, tx2 = tx2_id) %>%
    dplyr::select(tx1, tx2, everything())
  return(diff_df)
}

calculateAllPairwiseDiffs <- function(transcript_ids, exons_list){
  print(transcript_ids[1])
  
  #Extract selected transcripts
  tx_list = exons_list[transcript_ids]
  
  #Construct all pairs of transcripts
  tx_pairs = combn(transcript_ids, 2, simplify = FALSE)
  
  #Construct all pairs of transcripts
  diff_df = purrr::map_df(tx_pairs, ~calculatePairDiff(.[1], .[2], tx_list))
  return(diff_df)
}

#Import revised transcript annotations
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")

#Import gene metadata
revised_gene_metadata = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2()
group_map = dplyr::select(revised_gene_metadata, transcript_id, group_id)

#Import QTLs
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
acldl_df = readRDS("results/trQTLs/variance_explained/acldl_compiled_varExp.rds")
qtl_df = bind_rows(salmonella_df, acldl_df)

#Identify promoter events
qtls = dplyr::filter(qtl_df, quant == "reviseAnnotations") %>% 
  dplyr::rename(transcript_id = phenotype_id) %>%
  dplyr::left_join(group_map) %>%
  tidyr::separate(group_id, c("ensembl_gene_id", "grp", "position"), sep = "\\.", remove = FALSE) %>%
  dplyr::filter(position == "upstream") %>%
  dplyr::select(group_id) %>%
  dplyr::distinct()

#Extract all transcripts
selected_transcripts = dplyr::filter(group_map, group_id %in% qtls$group_id) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(transcript_count = length(transcript_id)) %>%
  dplyr::filter(transcript_count <= 10) 
all_diffs = purrrlyr::by_slice(selected_transcripts, 
        ~calculateAllPairwiseDiffs(.$transcript_id, revised_granges), .collate = "rows")

#Summarize at the group level
group_summary = dplyr::group_by(all_diffs, group_id) %>% 
  dplyr::summarize(upstream = sum(upstream), downstream = sum(downstream), contained = sum(contained))
saveRDS(group_summary, "results/reviseAnnotations/promoter_events.rds")


#Identify 3UTR events
qtls = dplyr::filter(qtl_df, quant == "reviseAnnotations") %>% 
  dplyr::rename(transcript_id = phenotype_id) %>%
  dplyr::left_join(group_map) %>%
  tidyr::separate(group_id, c("ensembl_gene_id", "grp", "position"), sep = "\\.", remove = FALSE) %>%
  dplyr::filter(position == "downstream") %>%
  dplyr::select(group_id) %>%
  dplyr::distinct()

#Extract all transcripts
selected_transcripts = dplyr::filter(group_map, group_id %in% qtls$group_id) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(transcript_count = length(transcript_id)) %>%
  dplyr::filter(transcript_count <= 10) 
all_diffs = purrrlyr::by_slice(selected_transcripts, 
                               ~calculateAllPairwiseDiffs(.$transcript_id, revised_granges), .collate = "rows")

#Summarize at the group level
group_summary = dplyr::group_by(all_diffs, group_id) %>% 
  dplyr::summarize(upstream = sum(upstream), downstream = sum(downstream), contained = sum(contained))
saveRDS(group_summary, "results/reviseAnnotations/UTR_events.rds")



##### Tests cases for new promoter code ####
#Identify all promoter events
promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "IRF5") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts)
plotTranscripts(new_transcripts)


promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "KIAA1671") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts)
plotTranscripts(new_transcripts)


promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "IRF2") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts)
plotTranscripts(new_transcripts)


promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "CD40") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts)
plotTranscripts(new_transcripts)


promoter_events = dplyr::filter(revised_gene_metadata, gene_name == "SLC9B2") %>% 
  dplyr::filter(group_id %like% "upstream") %>% 
  dplyr::filter(group_id %like% "grp_1")

transcripts = revised_granges[promoter_events$transcript_id]
plotTranscripts(transcripts)

new_transcripts = fillMissingInternalExons(transcripts)
plotTranscripts(new_transcripts)

