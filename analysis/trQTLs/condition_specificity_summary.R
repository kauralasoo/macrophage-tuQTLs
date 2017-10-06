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

#Extract gene names for all QTL types
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
gene_names = rowData(se_featureCounts) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = gene_id, gene_id)

se_revised = readRDS("results/SummarizedExperiments/acLDL_salmon_reviseAnnotations.rds")
revised_gene_names = rowData(se_revised) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

se_ensembl = readRDS("results/SummarizedExperiments/acLDL_salmon_Ensembl_87.rds")
ensembl_gene_names = rowData(se_ensembl) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id)

#Salmonella
leafcutter_gene_names = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds") %>%
  rowData() %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

#AcLDL
acldl_leafcutter_gene_names = readRDS("results/SummarizedExperiments/acLDL_leafcutter_counts.rds") %>%
  rowData() %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

#Merge all of them together
salmonella_gene_names = dplyr::bind_rows(gene_names, revised_gene_names, ensembl_gene_names, leafcutter_gene_names)
acldl_gene_names = dplyr::bind_rows(gene_names, revised_gene_names, ensembl_gene_names, acldl_leafcutter_gene_names)

#Caclulate fractions
salmonella_fractions = purrr::map(salmonella_varexp_list, ~purrr::map(.,~dplyr::mutate(.,interaction_fraction = interaction/(genotype+interaction))))
acLDL_fractions = purrr::map(acLDL_varexp_list, ~purrr::map(.,~dplyr::mutate(.,interaction_fraction = interaction/(genotype+interaction))))

#Count condition-specific factions
salmonella_df = purrr::map(salmonella_fractions, ~purrr::map_df(., identity, .id = "condition")) %>% 
  purrr::map_df(identity, .id = "quant") %>%
  dplyr::left_join(salmonella_gene_names, by = "phenotype_id")
saveRDS(salmonella_df, "results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")

acldl_df = purrr::map(acLDL_fractions, ~purrr::map_df(., identity, .id = "condition")) %>% 
  purrr::map_df(identity, .id = "quant") %>%
  dplyr::left_join(acldl_gene_names, by = "phenotype_id")
saveRDS(acldl_df, "results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")


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
se_tpm = readRDS("results/SummarizedExperiments/salmonella_salmon_gene_abundances.rds")
sample_meta = colData(se_tpm) %>% tbl_df2()
tpm_matrix = assays(se_tpm)$tpms

#Calculcate mean expression of each gene in each conditon
mean_tpm = calculateMean(tpm_matrix, design = as.data.frame(sample_meta), factor = "condition_name")
tpm_df = dplyr::mutate(mean_tpm, gene_id = rownames(mean_tpm)) %>%
  dplyr::select(gene_id, everything()) %>%
  tidyr::gather("condition", "mean_tpm", IFNg:SL1344)

#Look at mean expression for leafcutter QTLs
a = dplyr::filter(salmonella_df, quant == "leafcutter", condition == "IFNg_SL1344") %>% 
  dplyr::left_join(salmonella_gene_names, by = "phenotype_id") %>%
  dplyr::mutate(condition = "naive") %>%
  dplyr::left_join(tpm_df, by = c("gene_id", "condition")) %>%
  dplyr::left_join(leaf_mean_df, by = c("phenotype_id", "condition")) %>%
  dplyr::mutate(response_QTL = ifelse(interaction_fraction > 0.8, TRUE, FALSE)) %>%
  dplyr::filter(!is.na(response_QTL))

ggplot(a, aes(x = log(mean_tpm + .1,2), color = response_QTL)) + geom_density()

#Revised
a = dplyr::filter(salmonella_df, quant == "reviseAnnotations", condition == "IFNg_SL1344") %>% 
  dplyr::left_join(salmonella_gene_names, by = "phenotype_id") %>%
  dplyr::mutate(condition = "naive") %>%
  dplyr::left_join(tpm_df, by = c("gene_id", "condition")) %>%
  dplyr::mutate(response_QTL = ifelse(interaction_fraction > 0.8, TRUE, FALSE)) %>%
  dplyr::filter(!is.na(response_QTL))

ggplot(a, aes(x = log(mean_tpm + .1,2), color = response_QTL)) + geom_density()


#Ensembl_87
a = dplyr::filter(salmonella_df, quant == "Ensembl_87", condition == "IFNg_SL1344") %>% 
  dplyr::left_join(salmonella_gene_names, by = "phenotype_id") %>%
  dplyr::mutate(condition = "naive") %>%
  dplyr::left_join(tpm_df, by = c("gene_id", "condition")) %>%
  dplyr::mutate(response_QTL = ifelse(interaction_fraction > 0.8, TRUE, FALSE)) %>%
  dplyr::filter(!is.na(response_QTL))

ggplot(a, aes(x = log(mean_tpm + .1,2), color = response_QTL)) + geom_density()


#featureCounts
a = dplyr::filter(salmonella_df, quant == "featureCounts", condition == "IFNg_SL1344") %>% 
  dplyr::mutate(gene_id = phenotype_id) %>%
  dplyr::mutate(condition = "naive") %>%
  dplyr::left_join(tpm_df, by = c("gene_id", "condition")) %>%
  dplyr::mutate(response_QTL = ifelse(interaction_fraction > 0.8, TRUE, FALSE)) %>%
  dplyr::filter(!is.na(response_QTL))

ggplot(a, aes(x = log(mean_tpm + .1,2), color = response_QTL)) + geom_density()


