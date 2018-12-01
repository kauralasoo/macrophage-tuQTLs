library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
load_all("analysis/housekeeping/")
load_all("../seqUtils/")

#Import minimal QTL p-values
qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

#Extract QTLs
selected_qtls = qtls[c("leafcutter", "featureCounts", "Ensembl_87")]
naive_qtls = purrr::map(selected_qtls, ~.$naive)

#Combine txrevise QTLs
txrevise_qtls = dplyr::bind_rows(purrr::map_df(qtls$txrevise_promoters, identity, .id = "condition"),
                                 purrr::map_df(qtls$txrevise_contained, identity, .id = "condition"),
                                 purrr::map_df(qtls$txrevise_ends, identity, .id = "condition")) %>% 
  dplyr::filter(condition == "naive")
naive_qtls$reviseAnnotations = txrevise_qtls

#Add exptected p-values
qq_qtls = purrr::map_df(naive_qtls, ~dplyr::mutate(., p_eigen = p_beta) %>% addExpectedPvalue(.), .id = "quant") %>%
  dplyr::left_join(phenotypeFriendlyNames(), by = "quant")

#Count the number of phenotypes measured
pheno_count = dplyr::group_by(qq_qtls, phenotype) %>% dplyr::summarise(phenotype_count = length(phenotype_id))

qqplot = ggplot(qq_qtls, aes(x = -log(p_expected,10), y = -log(p_eigen, 10))) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "black") +
  facet_wrap(~phenotype) +
  xlab(expression(paste(Log[10], " expected p-value", sep = ""))) +
  ylab(expression(paste(Log[10], " permutation p-value (beta approximation)", sep = ""))) + 
  theme_light()

ggsave("results/figures/qtl_permutation_QQ_plots.png", plot = qqplot, width = 5, height = 5)
