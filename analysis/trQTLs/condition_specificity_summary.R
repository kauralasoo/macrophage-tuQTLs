library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import var explained for the Salmonella dataset
salmonella_varexp_list = list(Ensembl_87 = readRDS("results/trQTLs/variance_explained/salmonella_Ensembl_87_varExp.rds"),
                              reviseAnnotations = readRDS("results/trQTLs/variance_explained/salmonella_reviseAnnotations_varExp.rds"),
                              leafcutter = readRDS("results/trQTLs/variance_explained/salmonella_leafcutter_varExp.rds"),
                              featureCounts = readRDS("results/trQTLs/variance_explained/salmonella_featureCounts_varExp.rds"))
acLDL_varexp_list = list(Ensembl_87 = readRDS("results/trQTLs/variance_explained/acLDL_Ensembl_87_varExp.rds"),
                              reviseAnnotations = readRDS("results/trQTLs/variance_explained/acLDL_reviseAnnotations_varExp.rds"),
                              leafcutter = readRDS("results/trQTLs/variance_explained/acLDL_leafcutter_varExp.rds"),
                              featureCounts = readRDS("results/trQTLs/variance_explained/acLDL_featureCounts_varExp.rds"))                            

#Caclulate fractions
salmonella_fractions = purrr::map(salmonella_varexp_list, ~purrr::map(.,~dplyr::mutate(.,interaction_fraction = interaction/(genotype+interaction))))
acLDL_fractions = purrr::map(acLDL_varexp_list, ~purrr::map(.,~dplyr::mutate(.,interaction_fraction = interaction/(genotype+interaction))))

#Count condition-specific factions
salmonella_df = purrr::map(salmonella_fractions, ~purrr::map_df(., identity, .id = "condition")) %>% purrr::map_df(identity, .id = "quant") 
acldl_df = purrr::map(acLDL_fractions, ~purrr::map_df(., identity, .id = "condition")) %>% purrr::map_df(identity, .id = "quant") 

salmonella_fraction = dplyr::group_by(salmonella_df, quant, condition) %>% 
  dplyr::mutate(interaction_fraction = ifelse(is.na(interaction_fraction), 0, interaction_fraction)) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), interaction_count = sum(interaction_fraction > .5)) %>% 
  dplyr::mutate(interaction_percent = interaction_count/qtl_count)

acldl_fraction = dplyr::group_by(acldl_df, quant, condition) %>% 
  dplyr::mutate(interaction_fraction = ifelse(is.na(interaction_fraction), 0, interaction_fraction)) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), interaction_count = sum(interaction_fraction > .5)) %>% 
  dplyr::mutate(interaction_percent = interaction_count/qtl_count)

fraction_df = dplyr::bind_rows(salmonella_fraction, acldl_fraction) %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(phenotypeFriendlyNames())

#Make a barplot with the fraction of condition_specific QTLs
response_fraction_plot = ggplot(fraction_df, aes(x = phenotype, y = interaction_percent)) + 
  facet_wrap(~figure_name, nrow = 1) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Response QTL fraction") +
  xlab("Quantification method")
ggsave("results/figures/response_fraction_plot.pdf", plot = response_fraction_plot, width = 5, height = 3)


#Look at systematic differences in mean expression
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")

