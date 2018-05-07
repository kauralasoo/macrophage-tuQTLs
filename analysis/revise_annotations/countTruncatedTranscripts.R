library("dplyr")
library("devtools")
library("ggplot2")
load_all("../txrevise/")

#Import transcript metadata
tx_meta = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds") %>%
  tbl_df()

#Count truncated protein_coding transcripts
protein_tx = dplyr::filter(tx_meta, gene_biotype == "protein_coding", transcript_biotype == "protein_coding")
protein_count = dplyr::mutate(protein_tx, is_truncated = ifelse(cds_start_end_NF > 0, TRUE, FALSE)) %>% 
  dplyr::group_by(is_truncated) %>% 
  dplyr::summarise(transcript_count = length(ensembl_transcript_id)) %>%
  dplyr::mutate(type = "protein coding")

#Include processed and NMD transcripts as well
incomplete_transcripts = c("nonsense_mediated_decay","processed_transcript","retained_intron")
all_tx = dplyr::filter(tx_meta, gene_biotype == "protein_coding")
all_count = all_tx %>%
  dplyr::mutate(cds_start_NF = ifelse(transcript_biotype %in% incomplete_transcripts, 1, cds_start_NF), 
                cds_end_NF = ifelse(transcript_biotype %in% incomplete_transcripts, 1, cds_end_NF)) %>% 
  dplyr::mutate(cds_start_end_NF = pmax(cds_start_NF, cds_end_NF)) %>%
  dplyr::mutate(is_truncated = ifelse(cds_start_end_NF > 0, TRUE, FALSE)) %>% 
  dplyr::group_by(is_truncated) %>% 
  dplyr::summarise(transcript_count = length(ensembl_transcript_id)) %>%
  dplyr::mutate(type = "all")

merged_counts = dplyr::bind_rows(protein_count, all_count) %>%
  dplyr::mutate(type = factor(type, levels = c("protein coding", "all")))
truncated_fractions = ggplot(merged_counts, aes(x = type, y = transcript_count, fill = is_truncated)) + 
  geom_bar(stat = "identity") + 
  theme_light() + 
  xlab("transcript type") + 
  ylab("Number of transcripts")
ggsave("results/figures/truncated_transcript_fraction.pdf", plot = truncated_fractions, width = 4, height = 4)

#percentages
25314/(62325 + 25314)
95780/(62664 + 95780)



