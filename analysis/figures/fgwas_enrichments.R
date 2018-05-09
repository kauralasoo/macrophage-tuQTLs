library("dplyr")
library("readr")
library("tximport")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

file_list = list.files("processed/salmonella/fgwas/output/", recursive = T)

df = data_frame(file_path = file_list[grep("params", file_list)]) %>%
  tidyr::separate(file_path, c("model","quant", "cond_string"), sep ="\\/", remove = F) %>% 
  tidyr::separate(cond_string, c("condition_name","empty"), sep = ".params" ) %>%
  dplyr::select(-empty) %>%
  dplyr::mutate(full_path = file.path("processed/salmonella/fgwas/output", file_path))

#Make friendly names for the fgwas parameters
friendly_parameters = data_frame(
  parameter = c("promoter_ln",
                "fiveUTR_ln",
                "CDS_ln",
                "intron_ln",
                "threeUTR_ln",
                "polyA_site_ln",
                "atac_peak_ln",
                "eCLIP_splicing_ln"),
  parameter_name = factor(c("promoter","5' UTR", "coding", "intron","3' UTR", "poly(A)", "open\n chromatin", "splicing\n factor"),
                levels = c("promoter","5' UTR", "coding", "intron","3' UTR", "poly(A)", "open\n chromatin", "splicing\n factor")))


# Import files
file_paths = setNames(df$full_path, df$file_path)
all_data = purrr::map_df(file_paths, ~readr::read_delim(., delim = " ", col_types = "cddd"), .id = "file_path") %>%
  dplyr::mutate(CI_hi = ifelse(is.na(CI_hi), estimate + (estimate-CI_lo), CI_hi)) %>%
  dplyr::mutate(CI_lo = ifelse(is.na(CI_lo), estimate - (CI_hi-estimate), CI_lo)) %>%
  dplyr::left_join(df, by = "file_path") %>%
  left_join(conditionFriendlyNames()) %>%
  left_join(phenotypeFriendlyNames()) %>%
  dplyr::filter(!(parameter %in% c("eCLIP_both_ln","eCLIP_3end_ln"))) %>%
  dplyr::left_join(friendly_parameters)

### Compare different quant methods
joint_model_data = dplyr::filter(all_data, model == "atac_peak+polyA_site+fiveUTR+threeUTR+CDS+intron+promoter+eCLIP_3end+eCLIP_splicing+eCLIP_both")
data = dplyr::filter(joint_model_data, figure_name == "N", !is.na(phenotype)) %>%
  dplyr::filter(quant %in% c("Ensembl_87", "featureCounts","leafcutter","reviseAnnotations"))

methods_plot = ggplot(data, aes(x = estimate, y = phenotype, xmin = CI_lo, xmax = CI_hi)) + 
  geom_point() + facet_grid(parameter_name~.) +
  geom_errorbarh(aes(height = 0)) +
  theme_light() +
  theme(axis.title.y = element_blank()) +
  xlab("Log(enrichment)") +
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_cartesian(xlim = c(-.2,3.9)) +
  theme(panel.spacing = unit(0.1, "lines"))

ggsave("results/figures/fgwas_methods.pdf", plot = methods_plot, width = 3.2, height = 4)


#Compare different positioms
data = dplyr::filter(joint_model_data, figure_name == "N", quant %in% 
                       c("txrevise_contained","txrevise_ends","txrevise_promoters","leafcutter"))

methods_plot = ggplot(data, aes(x = estimate, y = phenotype, xmin = CI_lo, xmax = CI_hi)) + 
  geom_point() + facet_grid(parameter_name~.) +
  geom_errorbarh(aes(height = 0)) +
  theme_light() +
  theme(axis.title.y = element_blank()) +
  xlab("Log(enrichment)") +
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_cartesian(xlim = c(-.2,3.9)) +
  theme(panel.spacing = unit(0.1, "lines"))
ggsave("results/figures/fgwas_position.pdf", plot = methods_plot, width = 3.2, height = 4)



#Repeat the analysis for individual models

### Compare different quant methods
indiv_model_data = dplyr::filter(all_data, model != "atac_peak+polyA_site+fiveUTR+threeUTR+CDS+intron+promoter+eCLIP_3end+eCLIP_splicing+eCLIP_both")
data = dplyr::filter(indiv_model_data, figure_name == "I", !is.na(phenotype)) %>%
  dplyr::filter(quant %in% c("Ensembl_87", "featureCounts","leafcutter","reviseAnnotations"))

methods_plot = ggplot(data, aes(x = estimate, y = phenotype, xmin = CI_lo, xmax = CI_hi)) + 
  geom_point() + facet_grid(parameter_name~.) +
  geom_errorbarh(aes(height = 0)) +
  theme_light() +
  theme(axis.title.y = element_blank()) +
  xlab("Log(enrichment)") +
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_cartesian(xlim = c(-.2,4.5)) +
  theme(panel.spacing = unit(0.1, "lines"))


#Compare different positioms
data = dplyr::filter(indiv_model_data, figure_name == "I", quant %in% 
                       c("txrevise_contained","txrevise_ends","txrevise_promoters","leafcutter"))

methods_plot = ggplot(data, aes(x = estimate, y = phenotype, xmin = CI_lo, xmax = CI_hi)) + 
  geom_point() + facet_grid(parameter_name~.) +
  geom_errorbarh(aes(height = 0)) +
  theme_light() +
  theme(axis.title.y = element_blank()) +
  xlab("Log(enrichment)") +
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_cartesian(xlim = c(-.2,4.5)) +
  theme(panel.spacing = unit(0.1, "lines"))
ggsave("results/figures/fgwas_position.pdf", plot = methods_plot, width = 3.2, height = 4)

