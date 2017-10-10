library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")[c("Ensembl_87","reviseAnnotations","leafcutter","featureCounts")]
acLDL_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")

coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
salmonella_effects = dplyr::select(salmonella_df, quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

response_colocs = dplyr::left_join(coloc_df, salmonella_effects, by = c("quant", "phenotype_id", "snp_id")) %>% 
  dplyr::select(trait, quant, phenotype_id, snp_id, gene_name, interaction_fraction, p_fdr) %>% 
  dplyr::distinct()
View(dplyr::filter(response_colocs, is.na(interaction_fraction)))

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