library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Functions
countMultigeneQTLs <- function(pvalues_df, gene_name_map, fdr_threshold = 0.1, R2_threshold = 0.8){
  #Find shared QTLs
  hits = dplyr::filter(pvalues_df, p_fdr < fdr_threshold)
  field1 = dplyr::transmute(hits, gene_id = phenotype_id, snp_id)
  field2 = dplyr::transmute(hits, gene_id_2 = phenotype_id, snp_id, chr = snp_chr, pos = snp_start)
  revised_olaps = findGWASOverlaps(field1, field2, vcf_file, min_r2 = R2_threshold)
  
  #Count all unique genes with QTLs
  unique_genes = dplyr::left_join(hits, gene_name_map, by = "phenotype_id") %>% dplyr::select(gene_id) %>%
    unique()
  
  #Identfy multi-gene QTLs
  named_pairs = dplyr::filter(revised_olaps, gene_id != gene_id_2) %>% 
    dplyr::rename(phenotype_id = gene_id, phenotype_id_2 = gene_id_2) %>% 
    dplyr::left_join(gene_name_map, by = "phenotype_id") %>%
    dplyr::left_join(gene_name_map, by = c("phenotype_id_2" = "phenotype_id"))
  revised_multigene = dplyr::select(named_pairs,gene_id.x, gene_id.y, R2) %>% dplyr::filter(gene_id.x != gene_id.y)
  multigene_genes = unique(c(revised_multigene$gene_id.x, revised_multigene$gene_id.y))
  
  #Calculate ratio
  ratio = length(multigene_genes)/(nrow(unique_genes))
  
  return(list(R2 = revised_multigene, ratio = ratio))
}

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")

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
featureCounts_names = tbl_df2(rowData(se_featureCounts)) %>%
  dplyr::transmute(phenotype_id = gene_id, gene_name, gene_id)

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


#Count multi-gene QTLs per quantification method
fc_multi_list = purrr::map(salmonella_qtls$featureCounts[1:4], 
                           ~countMultigeneQTLs(.,featureCounts_names, fdr_threshold = 0.01, R2_threshold = 0.8))
lc_multi_list = purrr::map(salmonella_qtls$leafcutter[1:4], 
                           ~countMultigeneQTLs(.,leafcutter_names, fdr_threshold = 0.01, R2_threshold = 0.8))
re_multi_list = purrr::map(salmonella_qtls$reviseAnnotations[1:4], 
                           ~countMultigeneQTLs(.,revised_names, fdr_threshold = 0.01, R2_threshold = 0.8))
en_multi_list = purrr::map(salmonella_qtls$Ensembl_87[1:4], 
                           ~countMultigeneQTLs(.,ensembl_names, fdr_threshold = 0.01, R2_threshold = 0.8))

#Bind all estimates together
fc_df = purrr::map_df(fc_multi_list, ~data.frame(method = "featureCounts", ratio = .$ratio), identity)
lc_df = purrr::map_df(lc_multi_list, ~data.frame(method = "leafcutter", ratio = .$ratio), identity)
re_df = purrr::map_df(re_multi_list, ~data.frame(method = "reviseAnnotations", ratio = .$ratio), identity)
en_df = purrr::map_df(en_multi_list, ~data.frame(method = "Ensembl_87", ratio = .$ratio), identity)
ratio_df = dplyr::bind_rows(fc_df, lc_df, re_df, en_df)

multi_gene_qtl_plot = ggplot(ratio_df, aes(x = method, y = ratio)) + geom_boxplot() + geom_point()
ggsave("results/figures/multi-gene_qtl_fraction.pdf", plot = multi_gene_qtl_plot, width = 4, height = 4)



