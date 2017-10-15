library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")
coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
salmonella_effects = dplyr::select(salmonella_df, quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

#Link response QTLs to coloc results
response_colocs = dplyr::left_join(coloc_df, salmonella_effects, by = c("quant", "phenotype_id", "snp_id")) %>% 
  dplyr::select(trait, gwas_lead, quant, phenotype_id, snp_id, gene_name, interaction_fraction, p_fdr) %>% 
  dplyr::distinct()

#Identify unique genes
filtered_colocs = dplyr::filter(response_colocs, !is.na(p_fdr)) %>% 
  dplyr::group_by(trait, gwas_lead, quant, gene_name) %>% 
  dplyr::summarize(interaction_fraction = max(interaction_fraction), p_fdr = min(p_fdr)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(interaction_fraction)) %>%
  dplyr::mutate(is_response = ifelse(p_fdr < 0.1 & interaction_fraction > 0.5, TRUE, FALSE))

#What fraction of colocs show evidence for being a response QTL?
cond_summary = dplyr::group_by(filtered_colocs, quant) %>% 
  dplyr::summarize(coloc_count = length(gene_name), response_count = sum(is_response), response_fraction = response_count/coloc_count)


#### AcLDL ####
#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/acLDL_GWAS_coloc_hits.rds")
coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")
salmonella_effects = dplyr::select(salmonella_df, quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

#Link response QTLs to coloc results
response_colocs = dplyr::left_join(coloc_df, salmonella_effects, by = c("quant", "phenotype_id", "snp_id")) %>% 
  dplyr::select(trait, gwas_lead, quant, phenotype_id, snp_id, gene_name, interaction_fraction, p_fdr) %>% 
  dplyr::distinct()

#Identify unique genes
filtered_colocs = dplyr::filter(response_colocs, !is.na(p_fdr)) %>% 
  dplyr::group_by(trait, gwas_lead, quant, gene_name) %>% 
  dplyr::summarize(interaction_fraction = max(interaction_fraction), p_fdr = min(p_fdr)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(interaction_fraction)) %>%
  dplyr::mutate(is_response = ifelse(p_fdr < 0.1 & interaction_fraction > 0.5, TRUE, FALSE))

#What fraction of colocs show evidence for being a response QTL?
cond_summary = dplyr::group_by(filtered_colocs, quant) %>% 
  dplyr::summarize(coloc_count = length(gene_name), response_count = sum(is_response), response_fraction = response_count/coloc_count)
