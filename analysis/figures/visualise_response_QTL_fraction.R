library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Reimport trQTL data frames
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
acldl_df = readRDS("results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")

###### Estimate the fraction of QTLs that are condition specific #####
salmonella_fraction = dplyr::group_by(salmonella_df, quant, condition) %>% 
  dplyr::mutate(interaction_fraction = ifelse(is.na(interaction_fraction), 0, interaction_fraction)) %>% 
  dplyr::mutate(is_response_QTL = ifelse(interaction_fraction > .5 & p_fdr < 0.1, TRUE, FALSE)) %>%
  dplyr::summarise(qtl_count = length(phenotype_id), interaction_count = sum(is_response_QTL)) %>% 
  dplyr::mutate(interaction_percent = interaction_count/qtl_count)

acldl_fraction = dplyr::group_by(acldl_df, quant, condition) %>% 
  dplyr::mutate(interaction_fraction = ifelse(is.na(interaction_fraction), 0, interaction_fraction)) %>%
  dplyr::mutate(is_response_QTL = ifelse(interaction_fraction > .5 & p_fdr < 0.1, TRUE, FALSE)) %>%
  dplyr::summarise(qtl_count = length(phenotype_id), interaction_count = sum(is_response_QTL)) %>% 
  dplyr::mutate(interaction_percent = interaction_count/qtl_count)

fraction_df = dplyr::bind_rows(salmonella_fraction, acldl_fraction) %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(phenotypeFriendlyNames())

#Make a line plot
response_fraction_plot = ggplot(fraction_df, aes(x = phenotype, y = interaction_percent, color = figure_name, group = figure_name)) + 
  geom_point() + geom_line() +
  scale_color_manual(name = "condition", values = conditionPalette()) +
  theme_light() +
  ylab("Response QTL fraction") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1), axis.title.x = element_blank())
ggsave("results/figures/response_fraction_plot.pdf", plot = response_fraction_plot, width = 3.5, height = 2.5)

  
