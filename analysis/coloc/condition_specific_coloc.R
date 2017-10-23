library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("UpSetR")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")
coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
salmonella_effects = dplyr::select(salmonella_df, quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

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
  dplyr::mutate(is_response = ifelse(p_fdr < 0.1 & interaction_fraction > 0.2, TRUE, FALSE)) %>%
  dplyr::mutate(is_response = ifelse(is.na(is_response), FALSE, is_response)) 

#Identify response colocs
response_coloc_hits = dplyr::filter(filtered_colocs, is_response)

#### AcLDL ####
#Import coloc overlaps
salmonella_olaps = readRDS("results/coloc/acLDL_GWAS_coloc_hits.rds")
acldl_coloc_df = purrr::map_df(salmonella_olaps, identity, .id = "quant")

#Import QTL condition-specificity estimates
salmonella_df = readRDS("results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")
salmonella_effects = dplyr::select(salmonella_df, quant, condition, phenotype_id, snp_id, interaction_fraction, p_fdr)

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
  dplyr::mutate(is_response = ifelse(p_fdr < 0.1 & interaction_fraction > 0.2, TRUE, FALSE)) %>%
  dplyr::mutate(is_response = ifelse(is.na(is_response), FALSE, is_response)) 

acldl_response_coloc_hits = dplyr::filter(acldl_filtered_colocs, is_response)

#Extract summarized traits
summarised_traits = dplyr::bind_rows(coloc_df, acldl_coloc_df) %>%
  dplyr::select(trait, summarised_trait) %>% dplyr::distinct()

#Estimate the fraction of colocs that are condition specific
cond_specific_colocs = dplyr::bind_rows(filtered_colocs, acldl_filtered_colocs) %>%
  dplyr::left_join(summarised_traits) %>%
  dplyr::group_by(quant, summarised_trait, gene_name) %>%
  dplyr::mutate(has_response = as.logical(max(is_response))) %>%
  dplyr::select(quant, gene_name, summarised_trait, has_response) %>%
  dplyr::distinct() %>% dplyr::group_by(quant, has_response) %>% 
  dplyr::filter(gene_name != "FADS2") %>% #Remove FADS2 because it has too many associations
  dplyr::summarise(response_count = length(gene_name)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(phenotypeFriendlyNames())


cond_fraction = cond_specific_colocs %>%
  dplyr::mutate(has_response = ifelse(has_response, "response", "other")) %>%
  tidyr::spread(has_response, response_count) %>%
  dplyr::mutate(response_fraction = response/(other+response))


#Make a plot highilight cons-specific colocs
cond_speicifc_coloc_plot = ggplot(cond_specific_colocs, aes(x = phenotype, y = response_count, fill = has_response)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
  theme(legend.position = "top") + 
  ylab("Number of colocalisations") +
  labs(fill = "response QTL")
ggsave("results/figures/coloc_response_fraction.pdf",cond_speicifc_coloc_plot, width = 2.8, height = 3.8)




#Quantify coloc sharing between different phenotypes
all_colocs = dplyr::bind_rows(filtered_colocs, acldl_filtered_colocs) %>%
  dplyr::left_join(summarised_traits) %>%
  dplyr::left_join(phenotypeFriendlyNames()) %>%
  dplyr::filter(gene_name != "FADS2")

#Count overlaps by quantification strategy (by gene-trait pair)
unique_trait_gene_pairs = dplyr::select(all_colocs, phenotype, summarised_trait, gene_name) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "FCGR2A;RP11-25K21.6", "FCGR2A", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "VAMP8;VAMP5", "VAMP8", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "APOPT1;RP11-73M18.2", "APOPT1", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "FADS2;FEN1;FADS1", "FADS2", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "ST7L;WNT2B", "ST7L", gene_name)) %>%
  dplyr::mutate(gene_name = ifelse(gene_name == "FDPS;RUSC1-AS1", "FDPS", gene_name)) %>%
  dplyr::filter(!is.na(phenotype))

overlap_counts = dplyr::mutate(unique_trait_gene_pairs, id = gene_name) %>%
  tidyr::spread(phenotype, id) %>%
  dplyr::mutate_at(.vars = vars(3:6), 
                   .funs = function(x){ifelse(is.na(x), 0,1)}) %>%
  as.data.frame()

pdf("results/figures/coloc_GWAS_overlap_UpSetR.pdf", width = 4.5, height = 3.5, onefile = FALSE)
upset(as.data.frame(overlap_counts), sets = rev(c("read count", "transcript ratio", "Leafcutter", "txrevise")), 
      order.by = "freq", keep.order = TRUE)
dev.off()

