library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
load_all("analysis/housekeeping/")

#Import minimal QTL p-values
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
acldl_qtls = readRDS("results/trQTLs/acLDL_trQTL_min_pvalues.rds")

#Count QTLs in each condition
salmonella_qtl_count = purrr::map_df(salmonella_qtls, ~purrr::map_df(., ~dplyr::filter(., p_fdr < 0.1) %>% 
              dplyr::summarise(qtl_count = length(phenotype_id)), .id = "condition_name"), .id = "quant") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("naive","IFNg","SL1344","IFNg_SL1344", "NI","NS","NIS"))) %>%
  dplyr::mutate(quant = factor(quant, levels = c("featureCounts","tpm","leafcutter", "Ensembl_87" ,"reviseAnnotations", "reviseAnnotations_groupwise")))
acldl_qtl_count = purrr::map_df(acldl_qtls, ~purrr::map_df(., ~dplyr::filter(., p_fdr < 0.1) %>% 
              dplyr::summarise(qtl_count = length(phenotype_id)), .id = "condition_name"), .id = "quant") %>%
  dplyr::mutate(condition_name = factor(condition_name, levels = c("Ctrl","AcLDL","Diff"))) %>%
  dplyr::mutate(quant = factor(quant, levels = c("featureCounts","tpm","leafcutter", "Ensembl_87" ,"reviseAnnotations", "reviseAnnotations_groupwise")))

#Merge both counts
qtl_counts = dplyr::bind_rows(salmonella_qtl_count, acldl_qtl_count) %>% 
  dplyr::left_join(conditionFriendlyNames(), by = "condition_name") %>% 
  dplyr::filter(!is.na(figure_name)) %>%
  dplyr::left_join(phenotypeFriendlyNames(), by = "quant") %>%
  dplyr::filter(!is.na(phenotype))

#Make plots of the QTL counts
qtl_count_plot = ggplot(qtl_counts, aes(x = figure_name, y = qtl_count, fill = phenotype)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_light() + 
  xlab("Condition") + 
  ylab("Number of QTLs (FDR 10%)")
ggsave("results/figures/qtl_count_plot.pdf", plot = qtl_count_plot, width = 5.5, height = 3)
