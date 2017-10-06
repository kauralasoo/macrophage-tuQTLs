library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Helper functions
estimateConditionVarianceExplained <- function(qtl_df, trait_matrix, sample_metadata, 
                                               vcf_file, qtl_formula, interaction_fomula){
  #TODO: Add code to filete the vcf file to speed up interaction testing
  
  qtl_df = dplyr::select(qtl_df, phenotype_id, snp_id)
  res_df = purrrlyr::by_row(qtl_df, ~testInteractionLme4(.$phenotype_id, .$snp_id, trait_matrix, 
                                                         sample_metadata, vcf_file,
                                                         formula_qtl, formula_interaction, return_value = "all")$interaction_model %>%
                              varianceExplained(), .collate = "rows") %>%
    dplyr::select(-.row, -type)
  
  #Add proper name
  colnames(res_df)[colnames(res_df) == "condition_name:genotype"] = "interaction"
  return(res_df)
}

varExplainedWrapper <- function(qtl_df, trait_matrix, sample_meta, ...){
  estimateConditionVarianceExplained(qtl_df, trait_matrix, sample_meta, ...)
}

#Import all QTL calls
qtls = readRDS("results/trQTLs/acLDL_trQTL_min_pvalues.rds")

#Import genotypes
vcf_file = readRDS("results/genotypes/acLDL/imputed.70_samples.sorted.filtered.named.rds")

#Define formulas for interaction testing
formula_qtl = as.formula("expression ~ (1|genotype) + (1|condition_name)")
formula_interaction = as.formula("expression ~ (1|condition_name) + (1|condition_name:genotype) + (1|genotype)")

#Define condition pairs
condition_list = list(AcLDL = c("Ctrl", "AcLDL"))

##### Leafcutter #####
se_leafcutter = readRDS("results/SummarizedExperiments/acLDL_leafcutter_counts.rds")
leafcutter_by_cond = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_leafcutter))

#Extract feature matrices
leafcutter_mat_list = purrr::map(leafcutter_by_cond, ~assays(.)$tpm_ratios %>%
                                   replaceNAsWithRowMeans() %>%
                                   quantileNormaliseRows())
leafcutter_sample_meta = purrr::map(leafcutter_by_cond, ~colData(.) %>% tbl_df2())
leafcutter_gene_meta = rowData(se_leafcutter) %>% tbl_df2()
leafcutter_qtl_list = purrr::map(qtls$leafcutter[c("AcLDL")], ~dplyr::filter(.,p_fdr < 0.1))

#Estimate variance explained for all QTLs detected in stimulated conditons
leafcutter_var_explained = purrr::pmap(list(leafcutter_qtl_list, leafcutter_mat_list, leafcutter_sample_meta),
                varExplainedWrapper,
                vcf_file, qtl_formula, interaction_fomula)
saveRDS(leafcutter_var_explained, "results/trQTLs/variance_explained/acLDL_leafcutter_varExp.rds")


###### reviseAnnotations #####
se_reviseAnnotations = readRDS("results/SummarizedExperiments/acLDL_salmon_reviseAnnotations.rds")
revised_by_cond = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_reviseAnnotations))

#Extract feature matrices
revised_mat_list = purrr::map(revised_by_cond, ~assays(.)$tpm_ratios %>%
                                   replaceNAsWithRowMeans() %>%
                                   quantileNormaliseRows())
revised_sample_meta = purrr::map(revised_by_cond, ~colData(.) %>% tbl_df2())
revised_gene_meta = rowData(se_reviseAnnotations) %>% tbl_df2()
revised_qtl_list = purrr::map(qtls$reviseAnnotations[c("AcLDL")], ~dplyr::filter(.,p_fdr < 0.1))

#Estimate variance explained for all QTLs detected in stimulated conditons
revised_var_explained = purrr::pmap(list(revised_qtl_list, revised_mat_list, revised_sample_meta),
                                       varExplainedWrapper,
                                       vcf_file, qtl_formula, interaction_fomula)
saveRDS(revised_var_explained, "results/trQTLs/variance_explained/acLDL_reviseAnnotations_varExp.rds")


###### Ensembl_87 #####
se_ensembl = readRDS("results/SummarizedExperiments/acLDL_salmon_Ensembl_87.rds")
ensembl_by_cond = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_ensembl))

#Extract feature matrices
ensembl_mat_list = purrr::map(ensembl_by_cond, ~assays(.)$tpm_ratios %>%
                                replaceNAsWithRowMeans() %>%
                                quantileNormaliseRows())
ensembl_sample_meta = purrr::map(ensembl_by_cond, ~colData(.) %>% tbl_df2())
ensembl_gene_meta = rowData(se_ensembl) %>% tbl_df2()
ensembl_qtl_list = purrr::map(qtls$Ensembl_87[c("AcLDL")], ~dplyr::filter(.,p_fdr < 0.1))

#Estimate variance explained for all QTLs detected in stimulated conditons
ensembl_var_explained = purrr::pmap(list(ensembl_qtl_list, ensembl_mat_list, ensembl_sample_meta),
                                    varExplainedWrapper,
                                    vcf_file, qtl_formula, interaction_fomula)
saveRDS(ensembl_var_explained, "results/trQTLs/variance_explained/acLDL_Ensembl_87_varExp.rds")


###### featureCounts ######
se_featureCounts = readRDS("results/SummarizedExperiments/acLDL_featureCounts.rds")
fc_by_cond = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,se_featureCounts))

#Extract feature matrices
fc_mat_list = purrr::map(fc_by_cond, ~assays(.)$cqn)
fc_sample_meta = purrr::map(fc_by_cond, ~colData(.) %>% tbl_df2())
fc_gene_meta = rowData(se_featureCounts) %>% tbl_df2()
fc_qtl_list = purrr::map(qtls$featureCounts[c("AcLDL")], ~dplyr::filter(.,p_fdr < 0.1))

#Estimate variance explained for all QTLs detected in stimulated conditons
fc_var_explained = purrr::pmap(list(fc_qtl_list, fc_mat_list, fc_sample_meta),
                                       varExplainedWrapper,
                                       vcf_file, qtl_formula, interaction_fomula)
saveRDS(fc_var_explained, "results/trQTLs/variance_explained/acLDL_featureCounts_varExp.rds")





#Visualize a couple of examples
SPOPL_data = constructQtlPlotDataFrame("3:151219642:151261387:clu_8038", "rs1882017", 
                                       assays(leafcutter_by_cond$IFNg)$tpm_ratios, 
                                       vcf_file$genotypes, 
                                       leadfcutter_sample_meta$IFNg, 
                                       dplyr::mutate(leafcutter_gene_meta, gene_id = transcript_id)) %>%
  dplyr::mutate(genotype_text = genotype_value)
plotQtlCol(SPOPL_data)

#Visualize a couple of examples
SPOPL_data = constructQtlPlotDataFrame("ENSG00000109861", "rs11019479", 
                                       assays(se_featureCounts)$cqn, 
                                       vcf_file$genotypes, 
                                       fc_sample_meta$IFNg, 
                                       fc_gene_meta) %>%
  dplyr::mutate(genotype_text = genotype_value)
plotQtlCol(SPOPL_data)



