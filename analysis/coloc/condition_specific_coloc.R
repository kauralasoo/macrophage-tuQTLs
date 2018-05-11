library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("UpSetR")
library("data.table")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")
coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
contained_df = dplyr::filter(salmonella_df, quant == "reviseAnnotations", phenotype_id %like% 'contained') %>% 
  dplyr::mutate(quant = "txrevise_contained")
salmonella_effects = dplyr::bind_rows(salmonella_df, contained_df) %>% 
  dplyr::select(quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

#Link response QTLs to coloc results
response_colocs = dplyr::left_join(coloc_df, salmonella_effects, by = c("quant", "phenotype_id", "snp_id")) %>% 
  dplyr::select(trait, gwas_lead, quant, phenotype_id, snp_id, gene_name, interaction_fraction, p_fdr) %>% 
  dplyr::distinct()

#Identify unique genes
filtered_colocs = dplyr::filter(response_colocs, !is.na(p_fdr)) %>% 
  dplyr::group_by(trait, gwas_lead, quant, gene_name) %>% 
  dplyr::summarize(interaction_fraction = max(interaction_fraction), p_fdr = min(p_fdr)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(interaction_fraction)) %>%
  dplyr::mutate(is_response = ifelse(p_fdr < 0.1 & interaction_fraction > 0.5, TRUE, FALSE)) %>%
  dplyr::mutate(is_response = ifelse(is.na(is_response), FALSE, is_response)) 

#Quantify coloc sharing between different phenotypes
txrevise_combined = dplyr::filter(filtered_colocs, quant %in% c("txrevise_promoters", "txrevise_ends", "txrevise_contained")) %>%
  dplyr::mutate(quant = "reviseAnnotations")
filtered_colocs_txrevise = dplyr::bind_rows(filtered_colocs, txrevise_combined)

#Identify response colocs
response_coloc_hits = dplyr::filter(filtered_colocs_txrevise, is_response)

#### AcLDL ####
#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/acLDL_GWAS_coloc_hits.rds")
acldl_coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")
contained_df = dplyr::filter(salmonella_df, quant == "reviseAnnotations", phenotype_id %like% 'contained') %>% 
  dplyr::mutate(quant = "txrevise_contained")
salmonella_effects = dplyr::bind_rows(salmonella_df, contained_df) %>% 
  dplyr::select(quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

#Link response QTLs to coloc results
response_colocs = dplyr::left_join(acldl_coloc_df, salmonella_effects, by = c("quant", "phenotype_id", "snp_id")) %>% 
  dplyr::select(trait, gwas_lead, quant, phenotype_id, snp_id, gene_name, interaction_fraction, p_fdr) %>% 
  dplyr::distinct()

#Identify unique genes
acldl_filtered_colocs = dplyr::filter(response_colocs, !is.na(p_fdr)) %>% 
  dplyr::group_by(trait, gwas_lead, quant, gene_name) %>% 
  dplyr::summarize(interaction_fraction = max(interaction_fraction), p_fdr = min(p_fdr)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(interaction_fraction)) %>%
  dplyr::mutate(is_response = ifelse(p_fdr < 0.1 & interaction_fraction > 0.5, TRUE, FALSE)) %>%
  dplyr::mutate(is_response = ifelse(is.na(is_response), FALSE, is_response)) 

#Quantify coloc sharing between different phenotypes
txrevise_combined = dplyr::filter(acldl_filtered_colocs, quant %in% c("txrevise_promoters", "txrevise_ends", "txrevise_contained")) %>%
  dplyr::mutate(quant = "reviseAnnotations")
acldl_filtered_colocs_txrevise = dplyr::bind_rows(acldl_filtered_colocs, txrevise_combined)


acldl_response_coloc_hits = dplyr::filter(acldl_filtered_colocs_txrevise, is_response)

#Extract summarized traits
summarised_traits = dplyr::bind_rows(coloc_df, acldl_coloc_df) %>%
  dplyr::select(trait, summarised_trait) %>% dplyr::distinct()

#Estimate the fraction of colocs that are condition specific
cond_specific_colocs = dplyr::bind_rows(filtered_colocs_txrevise, acldl_filtered_colocs_txrevise) %>%
  dplyr::left_join(summarised_traits) %>%
  dplyr::group_by(quant, summarised_trait, gene_name) %>%
  dplyr::mutate(has_response = as.logical(max(is_response))) %>%
  dplyr::select(quant, gene_name, summarised_trait, has_response) %>%
  dplyr::distinct() %>% dplyr::group_by(quant, has_response) %>% 
  dplyr::filter(!(gene_name %in% c("FADS2","FADS2;FEN1;FADS1"))) %>% #Remove FADS2 because it has too many associations
  dplyr::summarise(response_count = length(gene_name)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(phenotypeFriendlyNames())

#Estimate the fraction of condition-specific colocalisations
cond_fraction = cond_specific_colocs %>%
  dplyr::mutate(has_response = ifelse(has_response, "response", "other")) %>%
  tidyr::spread(has_response, response_count) %>%
  dplyr::mutate(response = ifelse(is.na(response), 0, response)) %>%
  dplyr::mutate(response_fraction = response/(other+response))
write.table(cond_fraction, "results/tables/coloc_response_fraction.txt", sep = "\t", quote = F, row.names = F)

#Compare quant methods
coloc_response_plot = dplyr::filter(cond_fraction, quant %in% c("Ensembl_87","featureCounts","leafcutter","reviseAnnotations")) %>%
  ggplot(aes(x = phenotype, y = response_fraction)) + 
  geom_bar(stat = "identity") + 
  coord_flip(ylim = c(0,0.2)) +
  theme_light() +
  xlab("") +
  ylab("Fraction of colocalisations \n that are response QTLs")
ggsave("results/figures/coloc_response_fraction.pdf",coloc_response_plot, width = 3, height = 3)

#Compare positions
coloc_response_plot = dplyr::filter(cond_fraction, quant %in% c("leafcutter","txrevise_promoters","txrevise_contained","txrevise_ends")) %>%
  ggplot(aes(x = phenotype, y = response_fraction)) + 
  geom_bar(stat = "identity") + 
  coord_flip(ylim = c(0,0.3)) +
  theme_light() +
  xlab("") +
  ylab("Fraction of colocalisations \n that are response QTLs")
ggsave("results/figures/coloc_splicing_response_fraction.pdf",coloc_response_plot, width = 2.8, height = 3)



#Quantify coloc sharing between different phenotypes
all_colocs = dplyr::bind_rows(filtered_colocs_txrevise, acldl_filtered_colocs_txrevise) %>%
  dplyr::left_join(summarised_traits) %>%
  dplyr::left_join(phenotypeFriendlyNames())

#Count overlaps by quantification strategy (by gene-trait pair)
unique_trait_gene_pairs = dplyr::select(all_colocs, phenotype, summarised_trait, gene_name) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "FCGR2A;RP11-25K21.6", "FCGR2A", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "VAMP8;VAMP5", "VAMP8", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "APOPT1;RP11-73M18.2", "APOPT1", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "FADS2;FEN1;FADS1", "FADS2", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "ST7L;WNT2B", "ST7L", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "FDPS;RUSC1-AS1", "FDPS", gene_name)) %>%
  dplyr::filter(!is.na(phenotype)) %>%
  dplyr::filter(gene_name != "FADS2")

overlap_counts = dplyr::mutate(unique_trait_gene_pairs, id = gene_name) %>%
  tidyr::spread(phenotype, id) %>%
  dplyr::mutate_at(.vars = vars(3:9), 
                   .funs = function(x){ifelse(is.na(x), 0,1)}) %>%
  as.data.frame()

#Comapre different quantification methods
pdf("results/figures/coloc_GWAS_overlap_UpSetR.pdf", width = 3.5, height = 3, onefile = FALSE)
upset(as.data.frame(overlap_counts), sets = rev(c("read count", "transcript usage", "Leafcutter", "txrevise")), 
      order.by = "freq", keep.order = TRUE)
dev.off()

#Comapre only splicing methods
pdf("results/figures/coloc_GWAS_splicing_overlap_UpSetR.pdf", width = 3.5, height = 3, onefile = FALSE)
upset(as.data.frame(overlap_counts), sets = rev(c("Leafcutter", "promoters", "internal exons", "3' ends")), 
      order.by = "freq", keep.order = TRUE)
dev.off()

