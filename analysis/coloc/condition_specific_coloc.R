library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")
acLDL_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")

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