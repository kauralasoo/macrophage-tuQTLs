library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import eQTL summaries
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
acldl_df = readRDS("results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")
compiled_df = dplyr::bind_rows(salmonella_df, acldl_df)

#Estimate the fraction of trQTLs that are shared with eQTLs (R2 > 0.8)
sharing_df = dplyr::filter(compiled_df, quant != "featureCounts") %>% 
  dplyr::group_by(quant, condition) %>% 
  dplyr::mutate(is_shared = ifelse(R2 > .8, TRUE, FALSE)) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), shared_count = sum(is_shared)) %>% 
  dplyr::mutate(shared_fraction = shared_count/qtl_count) %>% 
  ungroup() %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(phenotypeFriendlyNames())

trQTL_eQTL_sharing = ggplot(sharing_df, aes(x = phenotype, y = shared_fraction)) + 
  geom_boxplot() + 
  geom_point() +
  scale_y_continuous(limits = c(0, .16)) + 
  ylab("Fraction of trQTLs shared with eQTLs (R2 > 0.8)") + 
  xlab("Quantification method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), axis.title.x = element_blank())
ggsave("results/figures/trQTL_eQTL_sharing.pdf", plot = trQTL_eQTL_sharing, width = 2, height = 4)


  