library("dplyr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Import event dataset
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
event_metadata = rowData(se_ensembl)
unique_genes = unique(event_metadata$ensembl_gene_id)

#Remove events on X and Y chromosomes
event_dataset = se_ensembl[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$transcript_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(event_dataset)) %>% 
  dplyr::filter(transcript_id %in% rownames(event_dataset)) %>%
  constructQTLtoolsGenePos()
output_path = "processed/salmonella/qtltools/input/Ensembl_87/"

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")

#Calculate fold-change
sample_list = colnames(event_conditions_renamed$naive)
fc_list = list(NI = assays(event_conditions_renamed$IFNg)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list],
               NS = assays(event_conditions_renamed$SL1344)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list],
               NIS = assays(event_conditions_renamed$IFNg_SL1344)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list])
normalised_fc_list = purrr::map(fc_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
qtltools_list = purrr::map(normalised_fc_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_list, output_path, file_suffix = "norm_prop")

#Calculate covariates for the fold-change
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_fc_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")





#Export reviseAnnotation quants for QTLtools
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
event_metadata = rowData(se_revised)
unique_genes = unique(event_metadata$gene_id)

#Remove events on X and Y chromosomes
event_dataset = se_revised[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$transcript_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(event_dataset)) %>% 
  dplyr::filter(transcript_id %in% rownames(event_dataset)) %>%
  constructQTLtoolsGenePos()
output_path = "processed/salmonella/qtltools/input/reviseAnnotations/"

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")

#Calculate fold-change
sample_list = colnames(event_conditions_renamed$naive)
fc_list = list(NI = assays(event_conditions_renamed$IFNg)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list],
               NS = assays(event_conditions_renamed$SL1344)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list],
               NIS = assays(event_conditions_renamed$IFNg_SL1344)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list])
normalised_fc_list = purrr::map(fc_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
qtltools_list = purrr::map(normalised_fc_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_list, output_path, file_suffix = "norm_prop")

#Calculate covariates for the fold-change
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_fc_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")





#Export LeafCutter counts for QTLtools
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
event_metadata = rowData(se_leafcutter)
unique_genes = unique(event_metadata$gene_id)

#Remove events on X and Y chromosomes
event_dataset = se_leafcutter[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$transcript_id,]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg","SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,event_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(event_dataset)) %>% 
  dplyr::filter(transcript_id %in% rownames(event_dataset)) %>%
  constructQTLtoolsGenePos()
output_path = "processed/salmonella/qtltools/input/leafcutter/"

#Extract ratio matrices
prop_list = purrr::map(event_conditions_renamed, ~assays(.)$tpm_ratios)

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
fastqtl_norm_prop_list = purrr::map(normalised_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(fastqtl_norm_prop_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")


#Calculate fold-change
sample_list = colnames(event_conditions_renamed$naive)
fc_list = list(NI = assays(event_conditions_renamed$IFNg)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list],
               NS = assays(event_conditions_renamed$SL1344)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list],
               NIS = assays(event_conditions_renamed$IFNg_SL1344)$tpm_ratios[,sample_list] -
                 assays(event_conditions_renamed$naive)$tpm_ratios[,sample_list])
normalised_fc_list = purrr::map(fc_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
qtltools_list = purrr::map(normalised_fc_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_list, output_path, file_suffix = "norm_prop")

#Calculate covariates for the fold-change
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_fc_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")





