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



#Export txrevise promoter quants for QTLtools
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_txrevise_promoters.rds")
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
output_path = "processed/salmonella/qtltools/input/txrevise_promoters/"

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





#Export reviseAnnotation quants for QTLtools (groupwise)
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
  dplyr::mutate(gene_id = group_id) %>%
  dplyr::filter(transcript_id %in% rownames(event_dataset)) %>%
  constructQTLtoolsGenePos()
output_path = "processed/salmonella/qtltools/input/reviseAnnotations_groupwise/"

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



#Filter reviseAnnotatons and export again
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

#Remove lowly expressed transcripts
event_conditions_filtered = purrr::map(event_conditions_renamed, removeLowlyExpressedTranscripts)

#Recalculate transcript ratios
event_metadata = purrr::map(event_conditions_filtered, ~rowData(.) %>% 
                              tbl_df2() %>%
                              dplyr::transmute(gene_id = group_id, transcript_id))
tpms = purrr::map(event_conditions_filtered, ~assays(.)$tpms)
prop_list = purrr::map2(tpms, event_metadata, ~calculateTranscriptRatios(.x,.y))

#Construct gene positions for QTL mapping
gene_pos_list = purrr::map(event_conditions_filtered, ~tbl_df2(rowData(.)) %>% 
  dplyr::filter(transcript_id %in% rownames(event_dataset)) %>%
  constructQTLtoolsGenePos())
output_path = "processed/salmonella/qtltools/input/reviseAnnotations_filtered/"

#Quantile normalise and remove NAs
normalised_list = purrr::map(prop_list, ~replaceNAsWithRowMeans(.) %>% quantileNormaliseRows())
fastqtl_norm_prop_list = purrr::map2(normalised_list, gene_pos_list, ~prepareQTLtoolsMatrix(.x,.y))
saveFastqtlMatrices(fastqtl_norm_prop_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(normalised_list, 
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



#Export Salmon TPMs for qtl mapping
se_gene = readRDS("results/SummarizedExperiments/salmonella_salmon_gene_abundances.rds")
event_metadata = rowData(se_gene)
unique_genes = unique(event_metadata$gene_id)

#Filtered gene metadata
filtered_transcscript_data = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.filtered.rds")

#Remove events on X and Y chromosomes
event_dataset = se_gene[event_metadata[!event_metadata$chr %in% c("X","Y","MT"),]$gene_id,]

#Filter by mean expression level
tpm_matrix = assays(event_dataset)$tpms
mean_by_condition = calculateMean(tpm_matrix, colData(event_dataset) %>% as.data.frame(), factor = "condition_name")
filtered_mean = mean_by_condition[rownames(mean_by_condition) %in% filtered_transcscript_data$ensembl_gene_id,]
expressed_genes = filtered_mean[apply(filtered_mean, 1, max) > 1,]
expressed_dataset = event_dataset[rownames(expressed_genes),]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg", "SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,expressed_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(expressed_dataset)) %>% 
  dplyr::filter(gene_id %in% rownames(expressed_dataset)) %>%
  dplyr::mutate(transcript_id = gene_id) %>%
  constructQTLtoolsGenePos()
output_path = "processed/salmonella/qtltools/input/tpm/"

#Extract ratio matrices
tpm_list = purrr::map(event_conditions_renamed, ~assays(.)$tpms)

#log transformation
log_tpm_list = purrr::map(tpm_list, ~log(. + 0.1,2))
qtltools_tpm_list = purrr::map(log_tpm_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_tpm_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(log_tpm_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")


#Calculate fold-change
sample_list = colnames(event_conditions_renamed$naive)
fc_list = list(NI = log_tpm_list$IFNg[,sample_list] -
                 log_tpm_list$naive[,sample_list],
               NS = log_tpm_list$SL1344[,sample_list] -
                 log_tpm_list$naive[,sample_list],
               NIS = log_tpm_list$IFNg_SL1344[,sample_list] -
                 log_tpm_list$naive[,sample_list])
qtltools_list = purrr::map(fc_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_list, output_path, file_suffix = "norm_prop")

#Calculate covariates for the fold-change
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(fc_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")


#Export featureCounts cqns for qtl mapping
se_gene = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
event_metadata = rowData(se_gene)
unique_genes = unique(event_metadata$gene_id)

#Identify expressed genes from Salmon TPMs
salmon_se = readRDS("results/SummarizedExperiments/salmonella_salmon_gene_abundances.rds")
salmon_meta = rowData(salmon_se)
salmon_dataset = salmon_se[salmon_meta[!salmon_meta$chr %in% c("X","Y","MT"),]$gene_id,]

tpm_matrix = assays(salmon_dataset)$tpms
mean_by_condition = calculateMean(tpm_matrix, colData(salmon_dataset) %>% as.data.frame(), factor = "condition_name")
expressed_genes = mean_by_condition[apply(mean_by_condition, 1, max) > 1,]

#Keep only expressed genes
expressed_dataset = se_gene[intersect(rownames(expressed_genes), unique_genes),]

#Extract lists for each condition
condition_list = idVectorToList(c("naive","IFNg", "SL1344","IFNg_SL1344"))
event_conditions = purrr::map(condition_list, ~extractConditionFromSummarizedExperiment(.,expressed_dataset))

#Rename columns
event_conditions_renamed = purrr::map(event_conditions, function(x){
  colnames(x) = x$genotype_id
  return(x)
})

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(expressed_dataset)) %>% 
  dplyr::filter(gene_id %in% rownames(expressed_dataset)) %>%
  dplyr::mutate(transcript_id = gene_id) %>%
  constructQTLtoolsGenePos()
output_path = "processed/salmonella/qtltools/input/featureCounts/"

#Extract cqn matrices
cqn_list = purrr::map(event_conditions_renamed, ~assays(.)$cqn)
qtltools_cqn_list = purrr::map(cqn_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_cqn_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(cqn_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")

#Calculate fold-change
sample_list = colnames(event_conditions_renamed$naive)
fc_list = list(NI = cqn_list$IFNg[,sample_list] -
                 cqn_list$naive[,sample_list],
               NS = cqn_list$SL1344[,sample_list] -
                 cqn_list$naive[,sample_list],
               NIS = cqn_list$IFNg_SL1344[,sample_list] -
                 cqn_list$naive[,sample_list])
qtltools_list = purrr::map(fc_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_list, output_path, file_suffix = "norm_prop")

#Calculate covariates for the fold-change
sample_meta = tbl_df2(colData(event_conditions_renamed$naive))
covariates_list = purrr::map(fc_list, 
                             ~performPCA(., sample_meta, n_pcs = 6, feature_id = "genotype_id")$pca_matrix %>%
                               dplyr::select(genotype_id, PC1, PC2, PC3, PC4, PC5, PC6) %>%
                               fastqtlMetadataToCovariates())
saveFastqtlMatrices(covariates_list, output_path, file_suffix = "covariates_prop")


