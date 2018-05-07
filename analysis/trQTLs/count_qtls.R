library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
load_all("analysis/housekeeping/")

#Import minimal QTL p-values
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
acldl_qtls = readRDS("results/trQTLs/acLDL_trQTL_min_pvalues.rds")

#Count QTLs for different parts of the gene
salmonella_qtl_count = purrr::map_df(salmonella_qtls, ~purrr::map_df(., ~dplyr::filter(., p_fdr < 0.1) %>% 
                                                                       dplyr::summarise(qtl_count = length(phenotype_id)), .id = "condition_name"), .id = "quant") %>%
  dplyr::filter(!(quant %in% c("reviseAnnotations_groupwise", "tpm")))
acldl_qtl_count = purrr::map_df(acldl_qtls, ~purrr::map_df(., ~dplyr::filter(., p_fdr < 0.1) %>% 
                                                             dplyr::summarise(qtl_count = length(phenotype_id)), .id = "condition_name"), .id = "quant") %>%
  dplyr::filter(!(quant %in% c("reviseAnnotations_groupwise", "tpm")))

#Merge both counts
qtl_counts = dplyr::bind_rows(salmonella_qtl_count, acldl_qtl_count) %>% 
  dplyr::left_join(conditionFriendlyNames(), by = "condition_name") %>% 
  dplyr::filter(!is.na(figure_name)) %>%
  dplyr::left_join(phenotypeFriendlyNames(), by = "quant") %>%
  dplyr::filter(!is.na(phenotype)) %>%
  dplyr::mutate(is_control = ifelse(figure_name %in% c("N","Ctrl"), TRUE, FALSE))

#Make plots of the QTL counts (main counts)
main_counts = dplyr::filter(qtl_counts, quant %in% c("Ensembl_87", "reviseAnnotations", "leafcutter", "featureCounts"))
qtl_count_plot = ggplot(main_counts, aes(x = phenotype, y = qtl_count, group = figure_name, color = figure_name, linetype = figure_name)) + 
  geom_point(position = position_dodge(0.05)) + 
  geom_line(position = position_dodge(0.05)) + 
  theme_light() + 
  coord_cartesian(ylim = c(0,3500)) +
  scale_color_manual(name = "condition", values = c("#636363","#67a9cf","#2166ac","#ef8a62","#bdbdbd","#b2182b")) +
  scale_linetype_manual(name = "condition", values = c("dotted", "solid", "solid", "solid","dotted", "solid")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
  ylab("Number of QTLs (FDR < 0.1)")

ggsave("results/figures/qtl_count_plot.pdf", plot = qtl_count_plot, width = 3.2, height = 2.8)



#Make plots of the QTL counts (tuQTL counts)
tuqtl_counts = dplyr::filter(qtl_counts, quant %in% c("leafcutter", "txrevise_promoters", "txrevise_contained", "txrevise_ends"))
qtl_count_plot = ggplot(tuqtl_counts, aes(x = phenotype, y = qtl_count, group = figure_name, color = figure_name, linetype = figure_name)) + 
  geom_point(position = position_dodge(0.05)) + 
  geom_line(position = position_dodge(0.05)) + 
  theme_light() + 
  coord_cartesian(ylim = c(0,2100)) +
  scale_color_manual(name = "condition", values = c("#636363","#67a9cf","#2166ac","#ef8a62","#bdbdbd","#b2182b")) +
  scale_linetype_manual(name = "condition", values = c("dotted", "solid", "solid", "solid","dotted", "solid")) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
  ylab("Number of QTLs (FDR < 0.1)")

ggsave("results/figures/qtl_count_plot_tuQTLs.pdf", plot = qtl_count_plot, width = 3.2, height = 2.8)




