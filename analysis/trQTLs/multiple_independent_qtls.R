library("dplyr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
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

countIndependentQTLs <- function(pvalues_list, gene_name_map, fdr_threshold = 0.1, R2_threshold = 0.2){
  
  #Extract QTLs
  qtls_list = purrr::map(pvalues_list, ~dplyr::filter(.,p_fdr < fdr_threshold) %>% 
                           dplyr::left_join(gene_name_map, by = "phenotype_id") %>%
                           dplyr::select(group_id, phenotype_id, gene_name, gene_id, p_nominal, snp_id) %>%
                           dplyr::arrange(gene_id, p_nominal) %>%
                           dplyr::group_by(gene_id) %>%
                           dplyr::mutate(qtl_count = length(gene_id)) %>%
                           dplyr::ungroup())
  
  #Count unique genes
  total_gene_count = purrr::map_df(qtls_list, identity, .id = "condition_name") %>%
    dplyr::select(gene_id, condition_name) %>%
    unique() %>%
    dplyr::group_by(condition_name) %>%
    dplyr::summarize(total_gene_count = length(gene_id))
  
  #Count multiple assoc for gene
  multi_gene_count = purrr::map_df(qtls_list, identity, .id = "condition_name") %>%
    dplyr::filter(qtl_count > 1) %>%
    dplyr::select(gene_id, condition_name) %>%
    unique() %>%
    dplyr::group_by(condition_name) %>%
    dplyr::summarize(multi_gene_count = length(gene_id))
  
  #Find indenpednent qtls
  independent_qtls = purrr::map(qtls_list, ~dplyr::select(., gene_id, snp_id, gene_name, phenotype_id) %>% 
                                  filterHitsR2(vcf_file$genotypes, R2_thresh = R2_threshold) %>% 
                                  dplyr::group_by(gene_id) %>%
                                  dplyr::mutate(qtl_count = length(gene_id)) %>%
                                  dplyr::ungroup() %>%
                                  dplyr::filter(qtl_count > 1) %>%
                                  dplyr::arrange(-qtl_count))
  
  independent_qtl_count = purrr::map_df(independent_qtls, identity, .id = "condition_name") %>%
    dplyr::select(gene_id, condition_name) %>%
    unique() %>%
    dplyr::group_by(condition_name) %>%
    dplyr::summarize(independent_gene_count = length(gene_id))
  
  multi_counts = dplyr::left_join(total_gene_count, multi_gene_count, by = "condition_name") %>%
    dplyr::left_join(independent_qtl_count, by = "condition_name") %>%
    dplyr::mutate(percent_all = independent_gene_count/total_gene_count) %>%
    dplyr::mutate(percent_multi = independent_gene_count/multi_gene_count)
  
  return(multi_counts)
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
revised_counts_02 = countIndependentQTLs(salmonella_qtls$reviseAnnotations[1:4], revised_names,
                                                  fdr_threshold = 0.1, R2_threshold = 0.2) %>%
  dplyr::mutate(R2 = 0.2)
revised_counts_08 = countIndependentQTLs(salmonella_qtls$reviseAnnotations[1:4], revised_names,
                                      fdr_threshold = 0.1, R2_threshold = 0.8) %>%
  dplyr::mutate(R2 = 0.8)
revised_counts = dplyr::bind_rows(revised_counts_02,revised_counts_08)
write.table(revised_counts, "results/figures/tables/revised_independent_qtls.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#### leafCutter ####
leafcutter_counts_02 = countIndependentQTLs(salmonella_qtls$leafcutter[1:4], leafcutter_names,
                                         fdr_threshold = 0.1, R2_threshold = 0.2) %>%
  dplyr::mutate(R2 = 0.2)
leafcutter_counts_08 = countIndependentQTLs(salmonella_qtls$leafcutter[1:4], leafcutter_names,
                                         fdr_threshold = 0.1, R2_threshold = 0.8) %>%
  dplyr::mutate(R2 = 0.8)
leafcutter_counts = dplyr::bind_rows(leafcutter_counts_02,leafcutter_counts_08)
write.table(leafcutter_counts, "results/figures/tables/leafcutter_independent_qtls.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)


##### Count QTLs that have multiple target genes ####

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



