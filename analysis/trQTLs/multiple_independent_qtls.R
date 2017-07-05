library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")

#Construct gene name map
revised_names = tbl_df2(rowData(se_revised)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id = ensembl_gene_id)
leafcutter_names = tbl_df2(rowData(se_leafcutter)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id = ensembl_gene_id)
ensembl_names = tbl_df2(rowData(se_ensembl)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id)

#Import QTLs
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

##### reviseAnnotations #####
#Identify genes with multiple independent QTLs
qtls = purrr::map_df(salmonella_qtls$reviseAnnotations[1:4], identity, .id = "condition_name") %>%
  dplyr::filter(p_fdr < 0.025) %>% 
  dplyr::left_join(revised_names, by = "phenotype_id") %>%
  dplyr::select(group_id, phenotype_id, gene_name, gene_id, p_nominal, snp_id) %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(qtl_count = length(gene_id)) %>%
  dplyr::ungroup() 

#Find indenpednent qtls
independent_qtls = dplyr::select(qtls, gene_id, snp_id, gene_name, phenotype_id) %>% 
  filterHitsR2(vcf_file$genotypes, R2_thresh = .2) %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(qtl_count = length(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(qtl_count > 1) %>%
  dplyr::arrange(-qtl_count)

a = length(unique(qtls$gene_id))
b = length(unique(independent_qtls$gene_id))
b/a

#### leafCutter ####
#Identify genes with multiple independent QTLs
qtls = purrr::map_df(salmonella_qtls$leafcutter[1:4], identity, .id = "condition_name") %>% 
  dplyr::filter(p_fdr < 0.025) %>% 
  dplyr::left_join(leafcutter_names, by = "phenotype_id") %>%
  dplyr::select(group_id, phenotype_id, gene_name, gene_id, p_nominal, snp_id) %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(qtl_count = length(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(gene_name))

#Find indenpednent qtls
independent_qtls = dplyr::select(qtls, gene_id, snp_id, gene_name, phenotype_id) %>% 
  filterHitsR2(vcf_file$genotypes, R2_thresh = .2) %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(qtl_count = length(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(qtl_count > 1) %>%
  dplyr::arrange(-qtl_count)

a = length(unique(qtls$gene_id))
b = length(unique(independent_qtls$gene_id))
b/a


#Ensembl
#### leafCutter ####
#Identify genes with multiple independent QTLs
qtls = purrr::map_df(salmonella_qtls$Ensembl_87[1:4], identity, .id = "condition_name") %>% 
  dplyr::filter(p_fdr < 0.025) %>% 
  dplyr::left_join(ensembl_names, by = "phenotype_id") %>%
  dplyr::select(group_id, phenotype_id, gene_name, gene_id, p_nominal, snp_id) %>%
  dplyr::arrange(gene_id, p_nominal) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(qtl_count = length(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(gene_name))

#Find indenpednent qtls
independent_qtls = dplyr::select(qtls, gene_id, snp_id, gene_name, phenotype_id) %>% 
  filterHitsR2(vcf_file$genotypes, R2_thresh = .2) %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(qtl_count = length(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(qtl_count > 1) %>%
  dplyr::arrange(-qtl_count)

a = length(unique(qtls$gene_id))
b = length(unique(independent_qtls$gene_id))
b/a


