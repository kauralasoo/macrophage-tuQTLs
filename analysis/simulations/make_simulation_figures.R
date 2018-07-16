library("dplyr")
library("devtools")
load_all("../seqUtils/")
library(ggplot2)

#Import transcript metadata
tx_meta = readRDS("results/simulations/transcript_meta.rds")
tx_diffs = readRDS("results/simulations/sim_one_diffs.rds")
all_tx = c(tx_diffs$full_tx, tx_diffs$truncated_tx)
gene_transcript_map = dplyr::transmute(tx_meta, gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id) %>%
  dplyr::filter(transcript_id %in% all_tx)

#extended_tx_diffs = readRDS("results/simulations/extended_transcript_diffs.rds") %>% 
#  dplyr::select(full_tx, truncated_tx, truncation)
#tx_diffs = readRDS("results/simulations/transcript_diffs.rds") %>% 
#  dplyr::select(-truncation) %>%
#  dplyr::left_join(extended_tx_diffs)

#Import QTLs
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

#Import lead variants
lead_variants = dplyr::filter(salmonella_qtls$Ensembl_87$naive, group_id %in% tx_meta$ensembl_gene_id) %>%
  dplyr::transmute(gene_id = group_id, snp_id)

#Import DEstatus
design = read.table("results/simulations/one_transcripts/sim_tx_info.txt", header = TRUE)
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
#sim_original = readRDS("processed/sim_original/matrices/original_transcripts.salmon_txrevise.rds")
sim_original = readRDS("processed/sim_both/matrices/both_transcripts.salmon_txrevise.rds")
tu_original = calculateTranscriptRatios(sim_original$abundance, gene_transcript_map)[selected_transcripts$transcript_id,]

#Import from extended transcripts
sim_extended = readRDS("processed/sim_one/matrices/both_transcripts.salmon_txrevise.rds")
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
upstream_original = readRDS("processed/sim_both/matrices/txrevise_upstream_both.salmon_txrevise.rds")
downstream_original = readRDS("processed/sim_both/matrices/txrevise_downstream_both.salmon_txrevise.rds")

#Construct event annotations
upstream_events = data_frame(transcript_id = rownames(upstream_original$abundance)) %>%
  tidyr::separate(transcript_id, c("gene_id", "suffix"), sep = "\\.grp_1", remove = FALSE) %>%
  dplyr::select(-suffix)

downstream_events = data_frame(transcript_id = rownames(downstream_original$abundance)) %>%
  tidyr::separate(transcript_id, c("gene_id", "suffix"), sep = "\\.grp_1", remove = FALSE) %>%
  dplyr::select(-suffix)

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
upstream_extended = readRDS("processed/sim_one/matrices/txrevise_upstream_both.salmon_txrevise.rds")
downstream_extended = readRDS("processed/sim_one/matrices/txrevise_downstream_both.salmon_txrevise.rds")

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
  dplyr::select(gene_id, upstream, downstream) %>%
  tidyr::gather(position, diff_length, upstream:downstream) %>% 
  dplyr::arrange(gene_id) %>%
  dplyr::mutate(is_truncated = ifelse(diff_length == 0, TRUE, FALSE)) %>%
  dplyr::mutate(event_type = ifelse(position == "upstream", "promoter", "3'end")) %>%
  dplyr::select(gene_id, diff_length, event_type, is_truncated)

extended_by_truncation = dplyr::left_join(gene_diffs, extended_effects, by = c("gene_id", "event_type"))

ggplot(extended_by_truncation, aes(x = estimate, fill = DEstatus)) + geom_histogram() + facet_grid(event_type~is_truncated)
ggplot(extended_by_truncation, aes(x = p.value, fill = DEstatus)) + geom_histogram() + facet_grid(event_type~is_truncated)

View(dplyr::filter(extended_by_truncation, event_type == "3'end", is_truncated == TRUE, DEstatus == TRUE))


plotTranscripts(exons[dplyr::filter(gene_transcript_map, gene_id == "ENSG00000033867")$transcript_id])
dplyr::filter(tx_diffs, truncated_tx == "ENST00000419036")

alt_events = readRDS("results/simulations/extended_tx_and_events.rds")
plotTranscripts(alt_events$ENSG00000033867$alt_events$ENSG00000033867.grp_1$upstream)
plotTranscripts(alt_events$ENSG00000033867$extended_tx$exons)


dplyr::filter(original_effects, transcript_id == "ENSG00000006283.grp_1.upstream.ENST00000514717")
dplyr::filter(extended_effects, transcript_id == "ENSG00000006283.grp_1.upstream.ENST00000514717")

dplyr::filter(original_effects, transcript_id == "ENSG00000163945.grp_1.upstream.ENST00000389851")
dplyr::filter(extended_effects, transcript_id == "ENSG00000163945.grp_1.upstream.ENST00000389851")


quant_alt_events = readRDS("results/simulations/qunatification_alt_events.rds")


ENSG00000175567

#Massive unexplained effect at the 3'UTR
dplyr::filter(original_effects, transcript_id == "ENSG00000183160.grp_1.downstream.ENST00000392806")
dplyr::filter(extended_effects, transcript_id == "ENSG00000183160.grp_1.downstream.ENST00000392806")

#Simulate data without alternatively spliced internal exons?
#Quantify without bias modelling?

#All constructed events
events = readRDS("results/simulations/one_both_alt_events.rds")

#Construct joint metadata
txrevise_metadata = rbind(upstream_events, downstream_events) %>% dplyr::arrange(gene_id) %>%
  dplyr::left_join(dplyr::transmute(tx_meta, gene_id = ensembl_gene_id, gene_name = external_gene_name) %>% dplyr::distinct()) %>% 
  dplyr::transmute(gene_id = transcript_id, gene_name)

#Make a QTL boxplot (relative expression)
selected_phenotype_id = "ENSG00000187147.grp_1.downstream.ENST00000440132"
selected_phenotype_id = "ENSG00000187147.grp_1.upstream.ENST00000484745"
selected_snp_id = "rs12125215"
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      downstream_extended$counts, vcf_file$genotypes, sample_metadata, txrevise_metadata)
ggplot(plot_data, aes(x = genotype_value, y = norm_exp)) + geom_point()




selected_phenotype_id = "ENSG00000126001.grp_1.downstream.ENST00000422671"
selected_snp_id = "rs3055798"
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      event_usage, vcf_file$genotypes, sample_metadata, txrevise_metadata)
ggplot(plot_data, aes(x = genotype_value, y = norm_exp)) + geom_boxplot() + geom_point()



#Option 1: Align to geneome first with hisat2, then use salmon (or RSEM) to quantify
#Option 2: Construct alternative event annotations that contain 


