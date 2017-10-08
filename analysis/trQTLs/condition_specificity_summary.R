library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

##### Import var explained for the Salmonella dataset #####
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


##### Import interaction test p-values ####
salmonella_test_list = list(Ensembl_87 = readRDS("results/trQTLs/variance_explained/salmonella_Ensembl_87_interaction_test.rds"),
                            reviseAnnotations = readRDS("results/trQTLs/variance_explained/salmonella_reviseAnnotations_interaction_test.rds"),
                            leafcutter = readRDS("results/trQTLs/variance_explained/salmonella_leafcutter_interaction_test.rds"),
                            featureCounts = readRDS("results/trQTLs/variance_explained/salmonella_featureCounts_interaction_test.rds"))
acLDL_test_list = list(Ensembl_87 = readRDS("results/trQTLs/variance_explained/acLDL_Ensembl_87_interaction_test.rds"),
                         reviseAnnotations = readRDS("results/trQTLs/variance_explained/acLDL_reviseAnnotations_interaction_test.rds"),
                         leafcutter = readRDS("results/trQTLs/variance_explained/acLDL_leafcutter_interaction_test.rds"),
                         featureCounts = readRDS("results/trQTLs/variance_explained/acLDL_featureCounts_interaction_test.rds"))

#Merge all test tesults together
salmonella_test_df = purrr::map(salmonella_test_list, ~purrr::map_df(., identity, .id = "condition")) %>% 
  purrr::map_df(identity, .id = "quant") 
acLDL_test_df = purrr::map(acLDL_test_list, ~purrr::map_df(., identity, .id = "condition")) %>% 
  purrr::map_df(identity, .id = "quant")

###### Extract gene names for all QTL types #####
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
gene_names = rowData(se_featureCounts) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = gene_id, gene_id)

#Salmonella
revised_gene_names = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

ensembl_gene_names = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds") %>%
  rowData(.) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id)

leafcutter_gene_names = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds") %>%
  rowData() %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

#AcLDL
acLDL_revised_gene_names = readRDS("results/SummarizedExperiments/acLDL_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

acLDL_ensembl_gene_names = readRDS("results/SummarizedExperiments/acLDL_salmon_Ensembl_87.rds") %>%
  rowData(.) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id)

acldl_leafcutter_gene_names = readRDS("results/SummarizedExperiments/acLDL_leafcutter_counts.rds") %>%
  rowData() %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

#Merge all of them together
salmonella_gene_names = dplyr::bind_rows(gene_names, revised_gene_names, ensembl_gene_names, leafcutter_gene_names)
acldl_gene_names = dplyr::bind_rows(gene_names, acLDL_revised_gene_names, acLDL_ensembl_gene_names, acldl_leafcutter_gene_names)

##### Import differential expression results #####
log2FC_list = list(IFNg = readr::read_tsv("results/DE/naive_vs_IFNg_DESeq2_fold_change.txt.gz", col_types = c("ccdddddd")),
                   SL1344 = readr::read_tsv("results/DE/naive_vs_Salmonella_DESeq2_fold_change.txt.gz", col_types = c("ccdddddd")),
                   IFNg_SL1344 = readr::read_tsv("results/DE/naive_vs_IFNg+Salmonella_DESeq2_fold_change.txt.gz", col_types = c("ccdddddd")),
                   AcLDL = readr::read_tsv("results/DE/Ctrl_vs_AcLDL_DESeq2_fold_change.txt", col_types = c("ccdddddd")))
log2FC_df = purrr::map_df(log2FC_list, identity, .id = "condition") %>%
  dplyr::select(condition, gene_id, gene_name, log2FoldChange)


##### Import mean gene expression in the naive condition
#Look at systematic differences in mean expression
se_tpm = readRDS("results/SummarizedExperiments/salmonella_salmon_gene_abundances.rds")
sample_meta = colData(se_tpm) %>% tbl_df2()
tpm_matrix = assays(se_tpm)$tpms

#Calculcate mean expression of each gene in each conditon
mean_tpm = calculateMean(tpm_matrix, design = as.data.frame(sample_meta), factor = "condition_name")
tpm_df = dplyr::mutate(mean_tpm, gene_id = rownames(mean_tpm)) %>%
  dplyr::select(gene_id, everything()) %>%
  tidyr::gather("condition", "mean_tpm", IFNg:SL1344)
naive_tpm = dplyr::filter(tpm_df, condition == "naive") %>%
  dplyr::transmute(gene_id, naive_tpm = mean_tpm) %>%
  tbl_df()

##### Merge all results together #####
salmonella_df = purrr::map(salmonella_fractions, ~purrr::map_df(., identity, .id = "condition")) %>% 
  purrr::map_df(identity, .id = "quant") %>%
  dplyr::left_join(salmonella_gene_names, by = "phenotype_id") %>%
  dplyr::left_join(salmonella_test_df, by = c("quant","condition", "phenotype_id","snp_id")) %>%
  dplyr::left_join(log2FC_df, by = c("condition", "gene_id")) %>%
  dplyr::left_join(naive_tpm, by = "gene_id")

acldl_df = purrr::map(acLDL_fractions, ~purrr::map_df(., identity, .id = "condition")) %>% 
  purrr::map_df(identity, .id = "quant") %>%
  dplyr::left_join(acldl_gene_names, by = "phenotype_id") %>%
  dplyr::left_join(acLDL_test_df, by = c("quant","condition", "phenotype_id","snp_id")) %>%
  dplyr::left_join(log2FC_df, by = c("condition", "gene_id")) %>%
  dplyr::left_join(naive_tpm, by = "gene_id")


##### Identify trQTLs that are linked to an eQTL of the same gene
#Salmonella
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")

#Extract eQTL lead variants
qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
salmonella_eQTLs = purrr::map_df(qtls$featureCounts, identity, .id = "condition") %>% 
  dplyr::transmute(gene_id = phenotype_id, eQTL_snp_id = snp_id, condition)

#Find linked eQTLs
trQTLs = dplyr::filter(salmonella_df, quant %in% c("reviseAnnotations", "Ensembl_87","leafcutter")) %>% 
  dplyr::transmute(condition, quant, phenotype_id, gene_id, snp_id) %>% 
  dplyr::left_join(salmonella_eQTLs, by = c("gene_id", "condition")) %>% 
  dplyr::filter(!is.na(eQTL_snp_id)) 
snp_pairs = dplyr::select(trQTLs, snp_id, eQTL_snp_id) %>% unique()
pair_r2 = purrrlyr::by_row(snp_pairs, ~calculatePairR2(.$snp_id, .$eQTL_snp_id, vcf_file$genotypes), .to = "R2", .collate = "rows")
salmonella_r2_df = dplyr::left_join(trQTLs, pair_r2, by = c("snp_id", "eQTL_snp_id")) %>% 
  dplyr::select(condition, quant, phenotype_id, snp_id, R2)

#AcLDL
vcf_file = readRDS("results/genotypes/acLDL/imputed.70_samples.sorted.filtered.named.rds")

qtls = readRDS("results/trQTLs/acLDL_trQTL_min_pvalues.rds")
acLDL_eQTLs = purrr::map_df(qtls$featureCounts, identity, .id = "condition") %>% 
  dplyr::transmute(gene_id = phenotype_id, eQTL_snp_id = snp_id, condition)

trQTLs = dplyr::filter(acldl_df, quant %in% c("reviseAnnotations", "Ensembl_87","leafcutter")) %>% 
  dplyr::transmute(condition, quant, phenotype_id, gene_id, snp_id) %>% 
  dplyr::left_join(acLDL_eQTLs, by = c("gene_id", "condition")) %>% 
  dplyr::filter(!is.na(eQTL_snp_id))
snp_pairs = dplyr::select(trQTLs, snp_id, eQTL_snp_id) %>% unique()
acLDL_pair_r2 = purrrlyr::by_row(snp_pairs, ~calculatePairR2(.$snp_id, .$eQTL_snp_id, vcf_file$genotypes), .to = "R2", .collate = "rows")
acLDL_r2_df = dplyr::left_join(trQTLs, acLDL_pair_r2, by = c("snp_id", "eQTL_snp_id")) %>% 
  dplyr::select(condition, quant, phenotype_id, snp_id, R2)

#Add R2 estimates to QTL df
salmonella_df_r2 = dplyr::left_join(salmonella_df, salmonella_r2_df, by = c("quant","condition","phenotype_id", "snp_id")) %>% 
  dplyr::mutate(R2 = ifelse(is.na(R2), 0, R2))
saveRDS(salmonella_df_r2, "results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")

acLDL_df_r2 = dplyr::left_join(acldl_df, acLDL_r2_df, by = c("quant","condition","phenotype_id", "snp_id")) %>% 
  dplyr::mutate(R2 = ifelse(is.na(R2), 0, R2))
saveRDS(acLDL_df_r2, "results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")

#Estimate the fraction of QTLs that are condition specific
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

#Make a barplot with the fraction of condition_specific QTLs
response_fraction_plot = ggplot(fraction_df, aes(x = phenotype, y = interaction_percent)) + 
  facet_wrap(~figure_name, nrow = 1) +
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Response QTL fraction") +
  xlab("Quantification method")
ggsave("results/figures/response_fraction_plot.pdf", plot = response_fraction_plot, width = 5, height = 3)


#Explore factors that explain response trQTLs
inter_genes = dplyr::filter(salmonella_df_r2, interaction_fraction > 0.5, p_fdr < 0.1, condition == "IFNg_SL1344")
ggplot(inter_genes, aes(x = log(naive_tpm + 1,2), color = quant)) + geom_histogram() + facet_wrap(~quant)
ggplot(inter_genes, aes(x = log(naive_tpm + 1,2), color = quant)) + geom_density()
ggplot(inter_genes, aes(x = log2FoldChange, color = quant)) + geom_histogram() + facet_wrap(~quant)
ggplot(inter_genes, aes(x = log2FoldChange, color = quant)) + geom_density(adjust = 2)




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


