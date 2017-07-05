library("dplyr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")

#Construct gene name map
revised_names = tbl_df2(rowData(se_revised)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id = ensembl_gene_id)
leafcutter_names = tbl_df2(rowData(se_leafcutter)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id = ensembl_gene_id)
ensembl_names = tbl_df2(rowData(se_ensembl)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id)

#Import QTLs
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

#Identify all genes with qtls
qtl_genes = dplyr::left_join(salmonella_qtls$reviseAnnotations$naive, revised_names, by = "phenotype_id") %>%
  dplyr::filter(p_fdr < 0.1) %>%
  dplyr::select(gene_id) %>%
  unique()

#naive QTLs
qtls = dplyr::left_join(salmonella_qtls$reviseAnnotations_groupwise$naive, revised_names, by = "phenotype_id") %>%
  tidyr::separate("group_id", c("gene1","group_name","position"), sep = "\\.", remove = FALSE) %>%
  dplyr::select(-gene1)

ifng_qtls = dplyr::left_join(salmonella_qtls$reviseAnnotations_groupwise$IFNg, revised_names, by = "phenotype_id") %>%
  tidyr::separate("group_id", c("gene1","group_name","position"), sep = "\\.", remove = FALSE) %>%
  dplyr::select(-gene1)

#Count QTLs per group
best_group = dplyr::group_by(qtls, gene_id, gene_name, group_name) %>% 
  dplyr::arrange(phenotype_id) %>%
  dplyr::filter(p_fdr < 0.1) %>% 
  dplyr::summarize(qtl_count = length(phenotype_id), min_p = min(p_beta)) %>% 
  ungroup() %>% dplyr::group_by(gene_id) %>% 
  arrange(gene_id, -qtl_count) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup()

best_filtered = dplyr::semi_join(qtls, best_group, by = c("gene_id", "group_name"))

#Perform Pi1 analysis
p_matrix = dplyr::select(best_filtered, gene_id, gene_name, position, p_beta) %>% 
  tidyr::spread(position, p_beta)

up_down = dplyr::filter(p_matrix, upstream < 1e-5, !is.na(downstream))
1 - qvalue::qvalue(up_down$downstream)$pi0

up_contained = dplyr::filter(p_matrix, upstream < 1e-5, !is.na(contained))
1 - qvalue::qvalue(up_contained$contained)$pi0

down_contained = dplyr::filter(p_matrix, downstream < 1e-5, !is.na(contained))
1 - qvalue::qvalue(down_contained$contained)$pi0

down_up = dplyr::filter(p_matrix, downstream < 1e-5, !is.na(upstream))
1 - qvalue::qvalue(down_up$upstream)$pi0

contained_up = dplyr::filter(p_matrix, contained < 1e-5, !is.na(upstream))
1 - qvalue::qvalue(contained_up$upstream)$pi0

contained_down = dplyr::filter(p_matrix, contained < 1e-5, !is.na(downstream))
1 - qvalue::qvalue(contained_down$downstream)$pi0


#Perform Pi1 between conditions
best_ifng = dplyr::semi_join(ifng_qtls, best_group, by = c("gene_id", "group_name"))

naive_hits = dplyr::filter(best_filtered, p_beta < 1e-5)
ifng_hits = dplyr::semi_join(ifng_qtls, naive_hits, by = "group_id")
1 - qvalue::qvalue(ifng_hits$p_beta)$pi0





