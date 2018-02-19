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
  dplyr::filter(!is.na(phenotype)) %>%
  dplyr::mutate(is_control = ifelse(figure_name %in% c("N","Ctrl"), TRUE, FALSE))

#Make plots of the QTL counts
qtl_count_plot = ggplot(qtl_counts, aes(x = phenotype, y = qtl_count, group = figure_name, color = figure_name, linetype = figure_name)) + 
  geom_point(position = position_dodge(0.05)) + 
  geom_line(position = position_dodge(0.05)) + 
  theme_light() + 
  coord_cartesian(ylim = c(0,3500)) +
  scale_color_manual(name = "condition", values = c("#636363","#67a9cf","#2166ac","#ef8a62","#bdbdbd","#b2182b")) +
  scale_linetype_manual(name = "condition", values = c("dotted", "solid", "solid", "solid","dotted", "solid")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
  ylab("Number of QTLs (FDR < 0.1)")

ggsave("results/figures/qtl_count_plot.pdf", plot = qtl_count_plot, width = 3.2, height = 2.8)
