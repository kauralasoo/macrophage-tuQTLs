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



#Perform tests for txrevise results


a = testInteraction(gene_id = "ENST00000514717", snp_id = "rs898455", tu_extended, sample_metadata, vcf_file, qtl_formula, interaction_formula, return_value = "all")



dplyr::filter(tx_meta, ensembl_transcript_id == "ENST00000514717")
original_transcripts = readRDS("processed/sim_extended/matrices/txrevise_downstream.salmon_txrevise.rds")
abundances_mat = original_transcripts$abundance
which(rownames(abundances_mat) %like% "ENSG00000006283")


