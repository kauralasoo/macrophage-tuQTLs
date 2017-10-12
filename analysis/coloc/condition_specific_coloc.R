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



#Partition into conditions
eqtl_coloc_counts = countConditionSpecificOverlaps(salmonella_olaps$tpm, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "RNA-seq")

#Partition into conditions
eqtl_coloc_counts = countConditionSpecificOverlaps(salmonella_olaps$leafcutter, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "RNA-seq")

#Partition into conditions
eqtl_coloc_counts = countConditionSpecificOverlaps(salmonella_olaps$revisedAnnotation, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "RNA-seq")

eqtl_coloc_counts = countConditionSpecificOverlaps(salmonella_olaps$ensembl_87, PP_power_thresh = 0.8, PP_coloc_thresh = .9)
eqtl_total_counts = group_by(eqtl_coloc_counts, figure_name) %>% 
  dplyr::summarise(overlap_count = sum(is_hit)) %>% 
  dplyr::mutate(total_overlap = cumsum(overlap_count)) %>%
  dplyr::mutate(phenotype = "RNA-seq")