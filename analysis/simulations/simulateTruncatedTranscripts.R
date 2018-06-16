library("BSgenome.Hsapiens.NCBI.GRCh38")
library("BSgenome")

#Import transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names=TRUE)
cdss = cdsBy(txdb, by = "tx", use.names=TRUE)

#Import QTL pairs
QTL_pairs = readRDS("results/simulations/trQTL_pair_diffs.rds")

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds"))
transcript_meta = dplyr::select(transcript_data, ensembl_transcript_id, cds_start_NF, cds_end_NF)

truncated_transcripts = dplyr::filter(truncated_transcripts, cds_start_NF == 1 | cds_end_NF == 1)

#Identify trQTL pairs with truncated transcripts
first_truncated = dplyr::semi_join(QTL_pairs, truncated_transcripts, by = c("tx1_id" = "ensembl_transcript_id"))
second_truncated = dplyr::semi_join(QTL_pairs, truncated_transcripts, by = c("tx2_id" = "ensembl_transcript_id"))
second_nonoverlap = dplyr::anti_join(second_truncated, first_truncated, by = "tx1_id")

#Make trancript pairs
truncated_pairs = dplyr::bind_rows(dplyr::transmute(first_truncated, full_tx = tx2_id, truncated_tx = tx1_id),
                 dplyr::transmute(second_nonoverlap, full_tx = tx1_id, truncated_tx = tx2_id)) %>%
  dplyr::left_join(truncated_transcripts, by = c("truncated_tx" = "ensembl_transcript_id")) %>%
  dplyr::mutate(truncation = case_when(
    cds_start_NF == 1 & cds_end_NF == 0 ~ "start",
    cds_start_NF == 0 & cds_end_NF == 1 ~ "end",
    cds_start_NF == 1 & cds_end_NF == 1 ~ "both"
  ))


#Calculate sequence differences in basepairs
findAllDiffs <- function(tx1, tx2, exons){
  print(paste(tx1, tx2))
  diff = txrevise::indentifyAddedRemovedRegions(tx1, tx2, exons) %>%
    calculateBasepairDifference()
}

#Find all differences between the two transcripts
tx1_list = as.list(truncated_pairs$full_tx)
tx2_list = as.list(truncated_pairs$truncated_tx)
all_differences = purrr::map2(tx1_list, tx2_list, ~findAllDiffs(.x, .y, exons)) %>% purrr::map_df(identity)

#Merge results
merged_diffs = dplyr::left_join(truncated_pairs, all_differences, by = c("full_tx" = "tx1_id")) %>% tbl_df()



tx_meta = dplyr::filter(transcript_data, ensembl_transcript_id %in% c("ENST00000514717", "ENST00000506520"))

gene_data = list()
gene_data$exons = exons[c("ENST00000514717", "ENST00000506520")]
gene_data$cdss = cdss[c("ENST00000514717", "ENST00000506520")]
gene_data$metadata = txrevise::filterTranscriptMetadata(tx_meta)

gene_extended_tx = txrevise::extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
gene_data_ext = txrevise::replaceExtendedTranscripts(gene_data, gene_extended_tx)

old_sequences = BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38, gene_data$exons)
new_sequences = BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38, gene_data_ext$exons)

old_fastas = DNAStringSet(lapply(old_sequences, unlist))
new_fastas = DNAStringSet(lapply(new_sequences, unlist))

wiggleplotr::plotTranscripts(gene_data$exons)
wiggleplotr::plotTranscripts(gene_data_ext$exons)






