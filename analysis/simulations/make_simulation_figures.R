library("dplyr")
library("devtools")
load_all("../seqUtils/")
library(ggplot2)

#Import transcript metadata
tx_meta = readRDS("results/simulations/transcript_meta.rds")
gene_transcript_map = dplyr::transmute(tx_meta, gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id)
tx_diffs = readRDS("results/simulations/transcript_diffs.rds")

#Import QTLs
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

#Import lead variants
lead_variants = dplyr::filter(salmonella_qtls$Ensembl_87$naive, group_id %in% tx_meta$ensembl_gene_id) %>%
  dplyr::transmute(gene_id = group_id, snp_id)

#Import DEstatus
design = read.table("results/simulations/original_transcripts/sim_tx_info.txt", header = TRUE)
de_status = dplyr::transmute(design, transcript_id = transcriptid, DEstatus = rowSums(design[,88:173])) %>% as_tibble() %>%
  dplyr::mutate(DEstatus = ifelse(DEstatus > 0, TRUE, FALSE))

#Extract selected transcripts
selected_transcripts = dplyr::group_by(gene_transcript_map, gene_id) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(lead_variants, by = "gene_id") %>%
  dplyr:::left_join(de_status, by = "transcript_id")

#Extract genotype
genotype_matrix = vcf_file$genotypes[lead_variants$snp_id,]

#Import salmon quantification results (original)
sim_original = readRDS("processed/sim_original/matrices/original_transcripts.salmon_txrevise.rds")
tu_original = calculateTranscriptRatios(sim_original$abundance, gene_transcript_map)[selected_transcripts$transcript_id,]

#Import from extended transcripts
sim_extended = readRDS("processed/sim_extended/matrices/original_transcripts.salmon_txrevise.rds")
tu_extended = calculateTranscriptRatios(sim_extended$abundance, gene_transcript_map)[selected_transcripts$transcript_id,]

#Construct sample metadata
sample_metadata = data_frame(sample_id = colnames(tu_original), genotype_id = colnames(genotype_matrix))

qtl_formula = as.formula(as.formula("expression ~ genotype"))
interaction_formula = as.formula(as.formula("expression ~ 1"))

#Perform all of the tests for original transcripts
snps_df = dplyr::transmute(selected_transcripts, phenotype_id = transcript_id, snp_id)
results = testMultipleInteractions(snps_df, tu_original, sample_metadata, vcf_file, qtl_formula, 
                                   interaction_formula, return_value = "all", id_field_separator = "-", lme4 = FALSE)
original_effects = purrr::map_df(results, ~broom::tidy(.$qtl_model), .id = "test_id") %>% 
  tidyr::separate(test_id, c("transcript_id", "snp_id"), sep = "-") %>% 
  dplyr::filter(term == "genotype") %>%
  dplyr::left_join(selected_transcripts, by = c("transcript_id", "snp_id")) %>%
  as_tibble() %>%
  dplyr::mutate(p_fdr = p.adjust(p.value, method = "fdr"))

#Make histograms
ggplot(original_effects, aes(x = p.value, fill = DEstatus)) + geom_histogram()
ggplot(original_effects, aes(x = estimate, fill = DEstatus)) + geom_histogram()


#Perform all of the tests for extended transripts
snps_df = dplyr::transmute(selected_transcripts, phenotype_id = transcript_id, snp_id)
results = testMultipleInteractions(snps_df, tu_extended, sample_metadata, vcf_file, qtl_formula, 
                                   interaction_formula, return_value = "all", id_field_separator = "-", lme4 = FALSE)
extended_effects = purrr::map_df(results, ~broom::tidy(.$qtl_model), .id = "test_id") %>% 
  tidyr::separate(test_id, c("transcript_id", "snp_id"), sep = "-") %>% 
  dplyr::filter(term == "genotype") %>%
  dplyr::left_join(selected_transcripts, by = c("transcript_id", "snp_id")) %>%
  as_tibble() %>%
  dplyr::mutate(p_fdr = p.adjust(p.value, method = "fdr"))

ggplot(extended_effects, aes(x = p.value, fill = DEstatus)) + geom_histogram()
ggplot(extended_effects, aes(x = estimate, fill = DEstatus)) + geom_histogram()



#### Perform tests for txrevise results
# Import quantification results
upstream_original = readRDS("processed/sim_original/matrices/txrevise_upstream.salmon_txrevise.rds")
downstream_original = readRDS("processed/sim_original/matrices/txrevise_downstream.salmon_txrevise.rds")

#Construct event annotations
upstream_events = data_frame(transcript_id = rownames(upstream_original$abundance)) %>%
  tidyr::separate(transcript_id, c("gene_id", "suffix"), sep = "\\.grp_1", remove = FALSE) %>%
  dplyr::select(-suffix)

downstream_events = data_frame(transcript_id = rownames(downstream_original$abundance)) %>%
  tidyr::separate(transcript_id, c("gene_id", "suffix"), sep = "\\.grp_1", remove = FALSE) %>%
  dplyr::select(-suffix) %>%
  dplyr::group_by(gene_id) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(event_type = "3'end")

selected_events = dplyr::bind_rows(dplyr::mutate(upstream_events, event_type = "promoter"), 
                                   dplyr::mutate(downstream_events, event_type = "3'end")) %>%
  dplyr::group_by(gene_id, event_type) %>% 
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(lead_variants, by = "gene_id") %>%
  dplyr:::left_join( dplyr::select(selected_transcripts, gene_id, DEstatus), by = "gene_id")

#Calculate event ratios
tu_up_original = calculateTranscriptRatios(upstream_original$abundance, upstream_events)
tu_down_original = calculateTranscriptRatios(downstream_original$abundance, downstream_events)
event_usage = rbind(tu_up_original, tu_down_original)[selected_events$transcript_id,]

#Run linear model
snps_df = dplyr::transmute(selected_events, phenotype_id = transcript_id, snp_id)
results = testMultipleInteractions(snps_df, event_usage, sample_metadata, vcf_file, qtl_formula, 
                                   interaction_formula, return_value = "all", id_field_separator = "-", lme4 = FALSE)
original_effects = purrr::map_df(results, ~broom::tidy(.$qtl_model), .id = "test_id") %>% 
  tidyr::separate(test_id, c("transcript_id", "snp_id"), sep = "-") %>% 
  dplyr::filter(term == "genotype") %>%
  dplyr::left_join(selected_events, by = c("transcript_id", "snp_id")) %>%
  as_tibble() %>%
  dplyr::mutate(p_fdr = p.adjust(p.value, method = "fdr")) %>% 
  dplyr::mutate(estimate = abs(estimate))

ggplot(original_effects, aes(x = p.value, fill = DEstatus)) + geom_histogram()
ggplot(original_effects, aes(x = estimate, fill = DEstatus)) + geom_histogram()


#repeat for extended transcripts
upstream_extended = readRDS("processed/sim_extended/matrices/txrevise_upstream.salmon_txrevise.rds")
downstream_extended = readRDS("processed/sim_extended/matrices/txrevise_downstream.salmon_txrevise.rds")

#Calculate event ratios
tu_up_original = calculateTranscriptRatios(upstream_extended$abundance, upstream_events)
tu_down_original = calculateTranscriptRatios(downstream_extended$abundance, downstream_events)
event_usage = rbind(tu_up_original, tu_down_original)[selected_events$transcript_id,]

#Run linear model
snps_df = dplyr::transmute(selected_events, phenotype_id = transcript_id, snp_id)
results = testMultipleInteractions(snps_df, event_usage, sample_metadata, vcf_file, qtl_formula, 
                                   interaction_formula, return_value = "all", id_field_separator = "-", lme4 = FALSE)
extended_effects = purrr::map_df(results, ~broom::tidy(.$qtl_model), .id = "test_id") %>% 
  tidyr::separate(test_id, c("transcript_id", "snp_id"), sep = "-") %>% 
  dplyr::filter(term == "genotype") %>%
  dplyr::left_join(selected_events, by = c("transcript_id", "snp_id")) %>%
  as_tibble() %>%
  dplyr::mutate(p_fdr = p.adjust(p.value, method = "fdr")) %>% 
  dplyr::mutate(estimate = abs(estimate))

ggplot(extended_effects, aes(x = p.value, fill = DEstatus)) + geom_histogram()
ggplot(extended_effects, aes(x = estimate, fill = DEstatus)) + geom_histogram()

#Identify truncated events
gene_diffs = dplyr::left_join(tx_diffs, gene_transcript_map, by = c("full_tx" = "transcript_id")) %>%
  dplyr::select(gene_id, upstream, downstream, truncation) %>%
  dplyr::filter(truncation != "both") %>%
  tidyr::gather(position, diff_length, upstream:downstream) %>%
  dplyr::mutate(event_type = ifelse(position == "upstream", "promoter", "3'end")) %>%
  dplyr::mutate(truncation = ifelse(truncation == "start", "promoter", "3'end")) %>%
  dplyr::mutate(is_truncated = ifelse(event_type == truncation, TRUE, FALSE)) %>%
  dplyr::filter(diff_length > 100) %>%
  dplyr::select(gene_id, diff_length, event_type, is_truncated)

extended_by_truncation = dplyr::left_join(gene_diffs, extended_effects, by = c("gene_id", "event_type"))

ggplot(extended_by_truncation, aes(x = estimate, fill = DEstatus)) + geom_histogram() + facet_grid(event_type~is_truncated)
ggplot(extended_by_truncation, aes(x = p.value, fill = DEstatus)) + geom_histogram() + facet_grid(event_type~is_truncated)

View(dplyr::filter(extended_by_truncation, event_type == "promoter", is_truncated == TRUE, DEstatus == TRUE))


plotTranscripts(exons[dplyr::filter(gene_transcript_map, gene_id == "ENSG00000033867")$transcript_id])
dplyr::filter(tx_diffs, truncated_tx == "ENST00000419036")

alt_events = readRDS("results/simulations/extended_tx_and_events.rds")
plotTranscripts(alt_events$ENSG00000033867$alt_events$ENSG00000033867.grp_1$upstream)
plotTranscripts(alt_events$ENSG00000033867$extended_tx$exons)


> dplyr::filter(original_effects, transcript_id == "ENSG00000006283.grp_1.upstream.ENST00000514717")
> dplyr::filter(extended_effects, transcript_id == "ENSG00000006283.grp_1.upstream.ENST00000514717")




