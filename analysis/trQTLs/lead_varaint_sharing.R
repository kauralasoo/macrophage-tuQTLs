library("devtools")
library("data.table")
library("dplyr")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

comparePairwise <- function(method_list, r2 = 0.8, fdr_thresh = 0.01){
  
  #Extract methods
  method_names = names(method_list)
  m1_table = method_list[[1]]
  m2_table = method_list[[2]]
  
  #Extract QTLs
  m1_qtls = dplyr::transmute(m1_table, gene_id, m1_snp = snp_id, m1_fdr = p_fdr) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(gene_id, m1_fdr) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup()
  
  m2_qtls = dplyr::transmute(m2_table, gene_id, m2_snp = snp_id, m2_fdr = p_fdr) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(gene_id, m2_fdr) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup()
  
  pairs = dplyr::left_join(m1_qtls, m2_qtls) %>% dplyr::filter(!(is.na(m2_snp))) %>% 
    dplyr::mutate(min_fdr = pmin(m1_fdr, m2_fdr)) %>%
    dplyr::filter(min_fdr < fdr_thresh)
  
  #Caclulate R2 for pairs
  pair_r2 = purrrlyr::by_row(pairs, ~calculatePairR2(.$m1_snp, .$m2_snp, vcf_file$genotypes), .to = "R2", .collate = "rows")
  
  #Estimate replication
  m1_table = table((dplyr::filter(pair_r2, m1_fdr < fdr_thresh) %>% dplyr::mutate(R2_thresh = R2 > r2))$R2_thresh)
  m1_rep = m1_table["TRUE"]/(m1_table["FALSE"]+m1_table["TRUE"])
  
  m2_table = table((dplyr::filter(pair_r2, m2_fdr < fdr_thresh) %>% dplyr::mutate(R2_thresh = R2 > r2))$R2_thresh)
  m2_rep = m2_table["TRUE"]/(m2_table["FALSE"]+m2_table["TRUE"])
  
  rep_df = data_frame(m1 = method_names, m2 = rev(method_names), replication = c(m1_rep, m2_rep))
}


#Sharing between Leafcutter QTLs and different txrevise QTLs
##### Identify trQTLs that are linked to an eQTL of the same gene

#Salmonella
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")

#Import gene names for all quant methods
leafcutter_gene_names = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds") %>%
  rowData() %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

revised_gene_names = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

ensembl_gene_names = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds") %>%
  rowData(.) %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id)


#Extract eQTL lead variants
qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
leafcutter_qtls = purrr::map_df(qtls$leafcutter, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>%
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr) %>%
  dplyr::left_join(leafcutter_gene_names)

txrevise_qtls = purrr::map_df(qtls$reviseAnnotations, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>%
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr) %>%
  dplyr::left_join(revised_gene_names)

ensembl_qtls = purrr::map_df(qtls$Ensembl_87, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>%
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr) %>%
  dplyr::left_join(ensembl_gene_names)

featureCounts_qtls = purrr::map_df(qtls$featureCounts, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>% 
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr, gene_id = group_id)


#Leafcutter vs txrevise
method_list = list(leafcutter = leafcutter_qtls, reviseAnnotations = txrevise_qtls)
leafcutter_vs_txrevise = comparePairwise(method_list)

#Leafcutter vs Ensembl 
method_list = list(leafcutter = leafcutter_qtls, Ensembl_87 = ensembl_qtls)
leafcutter_vs_ensembl = comparePairwise(method_list)

#Leafcutter vs featureCounts 
method_list = list(leafcutter = leafcutter_qtls, featureCounts = featureCounts_qtls)
leafcutter_vs_featureCounts = comparePairwise(method_list)

#Txrevise vs Ensembl 
method_list = list(reviseAnnotations = txrevise_qtls, Ensembl_87 = ensembl_qtls)
txrevise_vs_ensembl = comparePairwise(method_list)

#Txrevise vs featureCounts 
method_list = list(reviseAnnotations = txrevise_qtls, featureCounts = featureCounts_qtls)
txrevise_vs_featureCounts = comparePairwise(method_list)

#Ensembl vs featureCounts 
method_list = list(Ensembl_87 = ensembl_qtls, featureCounts = featureCounts_qtls)
ensmebl_vs_featureCounts = comparePairwise(method_list)

#Diagonal
diag_df = data_frame(m1 = unique(replication_df$m1), m2 = unique(replication_df$m1), replication = 1)

replication_df = dplyr::bind_rows(leafcutter_vs_txrevise,
                 leafcutter_vs_ensembl,
                 leafcutter_vs_featureCounts,
                 txrevise_vs_ensembl,
                 txrevise_vs_featureCounts,
                 ensmebl_vs_featureCounts, 
                 diag_df)

#Rename stuff
rep_df = dplyr::left_join(replication_df, phenotypeFriendlyNames(), by = c("m1" = "quant")) %>%
  dplyr::rename(p1 = phenotype) %>%
  dplyr::left_join(phenotypeFriendlyNames(), by = c("m2" = "quant")) %>%
  dplyr::rename(p2 = phenotype)

#Make a heatmap of sharing
qtl_replic_plot = ggplot(rep_df, aes(x = p2, y = p1, fill = replication, label = round(replication,2))) + geom_tile() +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", midpoint = 0, limits = c(0,1)) +
  theme_light() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  geom_text() +
  ylab("Query") +
  xlab("Replication") + 
  theme(legend.position = "none")
ggsave("results/figures/qtl_lead_variant_sharing.pdf", plot = qtl_replic_plot, width = 3, height = 2.5)



#Perform the same analysis for promoter vs UTR vs contained / Leafcutter
promoter_qtls = purrr::map_df(qtls$txrevise_promoters, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>%
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr) %>%
  dplyr::left_join(revised_gene_names)

end_qtls = purrr::map_df(qtls$txrevise_ends, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>%
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr) %>%
  dplyr::left_join(revised_gene_names)

middle_qtls = purrr::map_df(qtls$reviseAnnotations, identity, .id = "condition") %>% 
  dplyr::filter(condition == "naive") %>%
  dplyr::filter(phenotype_id %like% "contained") %>%
  dplyr::transmute(phenotype_id, snp_id, condition, p_fdr) %>%
  dplyr::left_join(revised_gene_names)

#Leafcutter vs Promoters 
method_list = list(leafcutter = leafcutter_qtls, txrevise_promoters = promoter_qtls)
leafcutter_vs_promoters = comparePairwise(method_list)

#Leafcutter vs Ends 
method_list = list(leafcutter = leafcutter_qtls, txrevise_ends = end_qtls)
leafcutter_vs_ends = comparePairwise(method_list)

#Leafcutter vs Middle 
method_list = list(leafcutter = leafcutter_qtls, txrevise_contained = middle_qtls)
leafcutter_vs_middle = comparePairwise(method_list)

event_df = dplyr::bind_rows(leafcutter_vs_promoters,leafcutter_vs_ends,leafcutter_vs_middle) %>%
  dplyr::filter(m2 == "leafcutter")

#Rename stuff
event_rename_df = dplyr::left_join(event_df, phenotypeFriendlyNames(), by = c("m1" = "quant")) %>%
  dplyr::rename(p1 = phenotype) %>%
  dplyr::left_join(phenotypeFriendlyNames(), by = c("m2" = "quant")) %>%
  dplyr::rename(p2 = phenotype)

event_replic_plot = ggplot(event_rename_df, aes(x = p2, y = p1, fill = replication, label = round(replication,2))) + geom_tile() +
  scale_fill_gradient2(space = "Lab", low = "#4575B4", mid = "#FFFFBF", high = "#E24C36", midpoint = 0, limits = c(0,1)) +
  theme_light() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text() +
  ylab("Query") +
  xlab("Replication") + 
  theme(legend.position = "none")
ggsave("results/figures/qtl_position_lead_variant_sharing.pdf", plot = event_replic_plot, width = 1.75, height = 2.5)



