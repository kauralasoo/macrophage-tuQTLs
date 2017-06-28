library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")

#Import minimal QTL p-values
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
acldl_qtls = readRDS("results/trQTLs/acLDL_trQTL_min_pvalues.rds")

#Count QTLs in each condition
salmonella_qtl_count = purrr::map_df(salmonella_qtls, ~purrr::map_df(., ~dplyr::filter(., p_fdr < 0.1) %>% 
              dplyr::summarise(qtl_count = length(phenotype_id)), .id = "condition"), .id = "quant") %>%
  dplyr::mutate(condition = factor(condition, levels = c("naive","IFNg","SL1344","IFNg_SL1344", "NI","NS","NIS"))) %>%
  dplyr::mutate(quant = factor(quant, levels = c("featureCounts","tpm","leafcutter", "Ensembl_87" ,"reviseAnnotations")))
acldl_qtl_count = purrr::map_df(acldl_qtls, ~purrr::map_df(., ~dplyr::filter(., p_fdr < 0.1) %>% 
              dplyr::summarise(qtl_count = length(phenotype_id)), .id = "condition"), .id = "quant") %>%
  dplyr::mutate(condition = factor(condition, levels = c("Ctrl","AcLDL","Diff"))) %>%
  dplyr::mutate(quant = factor(quant, levels = c("featureCounts","tpm","leafcutter", "Ensembl_87" ,"reviseAnnotations")))

#Make plots of the QTL counts
salmonella_plot = ggplot(salmonella_qtl_count, aes(x = condition, y = qtl_count, fill = quant)) + geom_bar(stat = "identity", position = "dodge")
ggsave("results/figures/salmonella_qtl_count.pdf", plot = salmonella_plot, width = 7, height = 4)

acLDL_plot = ggplot(acldl_qtl_count, aes(x = condition, y = qtl_count, fill = quant)) + geom_bar(stat = "identity", position = "dodge")
ggsave("results/figures/acLDL_qtl_count.pdf", plot = acLDL_plot, width = 5, height = 4)

