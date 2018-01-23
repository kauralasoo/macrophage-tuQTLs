library("devtools")
library("dplyr")
library("ggplot2")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import eQTL summaries
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
acldl_df = readRDS("results/trQTLs/variance_explained/acLDL_compiled_varExp.rds")
compiled_df = dplyr::bind_rows(salmonella_df, acldl_df)

#Estimate the fraction of trQTLs that are shared with eQTLs (R2 > 0.8)
sharing_df = dplyr::filter(compiled_df, quant != "featureCounts") %>% 
  dplyr::group_by(quant, condition) %>% 
  dplyr::mutate(is_shared = ifelse(R2 > .8, TRUE, FALSE)) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), shared_count = sum(is_shared)) %>% 
  dplyr::mutate(shared_fraction = shared_count/qtl_count) %>% 
  ungroup() %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(phenotypeFriendlyNames())

trQTL_eQTL_sharing = ggplot(sharing_df, aes(x = phenotype, y = shared_fraction)) + 
  geom_boxplot() + 
  geom_point() +
  scale_y_continuous(limits = c(0, .16)) + 
  ylab("Fraction of trQTLs shared with eQTLs (R2 > 0.8)") + 
  xlab("Quantification method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), axis.title.x = element_blank())
ggsave("results/figures/trQTL_eQTL_sharing.pdf", plot = trQTL_eQTL_sharing, width = 2, height = 4)


#Sharing between Leafcutter QTLs and different txrevise QTLs
##### Identify trQTLs that are linked to an eQTL of the same gene

#Salmonella
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")

leafcutter_gene_names = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds") %>%
  rowData() %>% tbl_df2() %>% 
  dplyr::transmute(phenotype_id = transcript_id, gene_id = ensembl_gene_id)

#Extract eQTL lead variants
qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
salmonella_leafcutter = purrr::map_df(qtls$leafcutter, identity, .id = "condition") %>% 
  dplyr::transmute(phenotype_id, eQTL_snp_id = snp_id, condition) %>%
  dplyr::left_join(leafcutter_gene_names) %>%
  dplyr::rename(lrasfcutter_id = phenotype_id)

#Find linked eQTLs
trQTLs = dplyr::filter(salmonella_df, quant %in% c("txrevise_promoters")) %>% 
  dplyr::transmute(condition, quant, phenotype_id, gene_id, snp_id) %>% 
  dplyr::left_join(salmonella_leafcutter, by = c("gene_id", "condition")) %>% 
  dplyr::filter(!is.na(eQTL_snp_id)) 
snp_pairs = dplyr::select(trQTLs, snp_id, eQTL_snp_id) %>% unique()
pair_r2 = purrrlyr::by_row(snp_pairs, ~calculatePairR2(.$snp_id, .$eQTL_snp_id, vcf_file$genotypes), .to = "R2", .collate = "rows")
salmonella_r2_df = dplyr::left_join(trQTLs, pair_r2, by = c("snp_id", "eQTL_snp_id")) %>% 
  dplyr::select(condition, quant, phenotype_id, snp_id, R2)
puqtl_counts = dplyr::group_by(salmonella_r2_df, phenotype_id) %>% dplyr::summarize(R2 = max(R2))
table(puqtl_counts$R2 > 0.8)
298/(298+819)

#Find linked eQTLs
sQTLs = dplyr::filter(salmonella_df, quant %in% c("reviseAnnotations"), phenotype_id %like% "contained") %>% 
  dplyr::transmute(condition, quant, phenotype_id, gene_id, snp_id) %>% 
  dplyr::left_join(salmonella_leafcutter, by = c("gene_id", "condition")) %>% 
  dplyr::filter(!is.na(eQTL_snp_id)) 
snp_pairs = dplyr::select(sQTLs, snp_id, eQTL_snp_id) %>% unique()
pair_r2 = purrrlyr::by_row(snp_pairs, ~calculatePairR2(.$snp_id, .$eQTL_snp_id, vcf_file$genotypes), .to = "R2", .collate = "rows")
contained_leafcutter_olap = dplyr::left_join(sQTLs, pair_r2, by = c("snp_id", "eQTL_snp_id")) %>% 
  dplyr::select(condition, quant, phenotype_id, snp_id, R2)
sqtl_counts = dplyr::group_by(contained_leafcutter_olap, phenotype_id) %>% dplyr::summarize(R2 = max(R2))
table(sqtl_counts$R2 > 0.8)
897/(1125+897)

apaQTLs = dplyr::filter(salmonella_df, quant %in% c("txrevise_ends")) %>% 
  dplyr::transmute(condition, quant, phenotype_id, gene_id, snp_id) %>% 
  dplyr::left_join(salmonella_leafcutter, by = c("gene_id", "condition")) %>% 
  dplyr::filter(!is.na(eQTL_snp_id)) 
snp_pairs = dplyr::select(apaQTLs, snp_id, eQTL_snp_id) %>% unique()
pair_r2 = purrrlyr::by_row(snp_pairs, ~calculatePairR2(.$snp_id, .$eQTL_snp_id, vcf_file$genotypes), .to = "R2", .collate = "rows")
apa_r2_df = dplyr::left_join(apaQTLs, pair_r2, by = c("snp_id", "eQTL_snp_id")) %>% 
  dplyr::select(condition, quant, phenotype_id, snp_id, R2)
apaqtl_counts = dplyr::group_by(apa_r2_df, phenotype_id) %>% dplyr::summarize(R2 = max(R2))
table(apaqtl_counts$R2 > 0.8)
274/(274+1101)


