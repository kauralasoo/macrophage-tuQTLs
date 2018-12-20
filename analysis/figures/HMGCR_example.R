library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("../wiggleplotr")
load_all("analysis/housekeeping/")
load_all("../txrevise/")

#Import SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")

#Extract TPM matrix
ratio_matrix = assays(se_revised)$tpm_ratios
tpm_matrix = assays(se_revised)$tpms
sample_meta = colData(se_revised) %>% tbl_df2()
gene_meta = rowData(se_revised) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = transcript_id)
gene_meta_ensembl = rowData(se_ensembl) %>% tbl_df2() 

#Import transcript coordinates
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names=TRUE)
cdss = cdsBy(txdb, by = "tx", use.names=TRUE)

#Import coloc resuts
coloc_hits = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
GRCh38_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
GRCh37_information = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")

#Extract HMGCR phenotype
selected_phenotype_id = "ENSG00000113161.grp_2.contained.ENST00000343975"
selected_snp_id = "rs3846662"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N"))

boxplot = plotQtlCol(plot_data) + ylab("txrevse: ENST00000343975 usage")
ggsave("results/figures/HMGCR_boxplot.pdf", plot = boxplot, width = 2.5, height = 3)


#Make a boxplot of the other trancript
selected_phenotype_id2 = "ENSG00000113161.grp_2.contained.ENST00000511206"

plot_data = constructQtlPlotDataFrame(selected_phenotype_id2, selected_snp_id, 
                                      ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N"))

boxplot2 = plotQtlCol(plot_data) + ylab("txrevise: ENST00000511206 usage")
ggsave("results/figures/HMGCR_boxplot2.pdf", plot = boxplot2, width = 2.5, height = 3)

#Make a read coverage plot
txrevise_ids = dplyr::filter(gene_meta, group_id == "ENSG00000113161.grp_2.contained")$transcript_id
ensembl_ids = filter(gene_meta_ensembl, gene_id == "ENSG00000113161")$transcript_id

selected_ensembl = c("ENST00000343975", "ENST00000287936")
selected_txrevise = c("ENSG00000113161.grp_2.contained.ENST00000511206", "ENSG00000113161.grp_2.contained.ENST00000343975")
plotTranscripts(revised_granges[txrevise_ids])
plotTranscripts(exons[selected_ensembl])


#Make read coverage plot as well
#Set up bigwig metadata
rna_meta_str1_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                                tbl_df2(colData(se_featureCounts)),
                                                "/Volumes/Ajamasin/bigwig/salmonella", 
                                                bigWig_suffix = ".str1.bw")
rna_meta_str2_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                                tbl_df2(colData(se_featureCounts)),
                                                "/Volumes/Ajamasin/bigwig/salmonella", 
                                                bigWig_suffix = ".str2.bw")


#Import full transcripts
ensembl_meta = rowData(se_ensembl) %>% tbl_df2() %>% 
  dplyr::filter(gene_name == "HMGCR") %>%
  dplyr::filter(transcript_biotype == "protein_coding") %>%
  dplyr::filter(transcript_id %in% c("ENST00000343975", "ENST00000287936"))
ensembl_exons = exons[ensembl_meta$transcript_id] %>% txrevise::removeMetadata()
ensembl_cdss = cdss[ensembl_meta$transcript_id] %>% txrevise::removeMetadata()
names(ensembl_exons) = paste("HMGCR: ", names(ensembl_exons), sep = "")
names(ensembl_cdss) = paste("HMGCR: ", names(ensembl_cdss), sep = "")

tx_ranges = revised_granges[selected_txrevise]
names(tx_ranges) = c("txrevise: ENST00000511206", "txrevise: ENST00000343975")


#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, selected_snp_id, vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(track_id = figure_name)

coverage_plot = plotCoverage(exons = c(ensembl_exons, rev(tx_ranges)), cdss = ensembl_cdss, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.5,0.5), 
                             return_subplots_list = F, coverage_type = "line", transcript_label = TRUE)
ggsave("results/figures/HMGCR_coverage_plot.pdf", plot = coverage_plot, width = 6, height = 4)




#Make a manhattan plot
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name", "type"))


##### Promoter QTL ####
qtl_df = data_frame(phenotype_id = selected_phenotype_id, snp_id = selected_snp_id, trait = "LDL")
qtl_paths = list(naive = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/reviseAnnotations/sorted/naive.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, 
                                         gwas_dir = "/Users/alasoo/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = GRCh37_information, GRCh38_variants = GRCh38_information, 
                                         cis_dist = 2e5, type = "QTLTools") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>% 
  dplyr::mutate(track_id = condition_name)


#Make a manhattan plot
region_coords = c(75250000, 75450000)
gwas_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id == "LDL"), region_coords, color_R2 = TRUE)
eqtl_data = dplyr::filter(qtl_summary, track_id != "LDL") %>% 
  dplyr::left_join(conditionFriendlyNames()) %>% 
  dplyr::mutate(track_id = figure_name)
eqtl_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(eqtl_data, track_id != "LDL"), region_coords, color_R2 = TRUE)


#Plot one transcript
tx_plot = plotTranscripts(exons = ensembl_exons[2], cdss = ensembl_cdss[2], rescale_introns = FALSE, region_coords = region_coords)


joint_plot = cowplot::plot_grid(gwas_manhattan, eqtl_manhattan, tx_plot,
                                align = "v", ncol = 1, rel_heights = c(3,3,3))
ggsave("results/figures/HMGCR_coloc.pdf", plot = joint_plot, width = 5, height = 4)




#Make boxplots of full-length transcript expression
#Extract TPM matrix
ensembl_ratio_matrix = assays(se_ensembl)$tpm_ratios
ensmebl_sample_meta = colData(se_ensembl) %>% tbl_df2()
ensembl_gene_meta = rowData(se_ensembl) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = transcript_id)

#Extract HMGCR phenotype
selected_phenotype_id = "ENST00000343975"
selected_snp_id = "rs3846662"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      ensembl_ratio_matrix, vcf_file$genotypes, ensmebl_sample_meta, ensembl_gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N"))

boxplot = plotQtlCol(plot_data) + ylab("ENST00000343975 usage") + scale_y_continuous(limits = c(0,0.5))
ggsave("results/figures/HMGCR_ENST00000343975_ensembl.pdf", plot = boxplot, width = 2.5, height = 3)

#Extract HMGCR phenotype
selected_phenotype_id = "ENST00000287936"
selected_snp_id = "rs3846662"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      ensembl_ratio_matrix, vcf_file$genotypes, ensmebl_sample_meta, ensembl_gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N"))

boxplot = plotQtlCol(plot_data) + ylab("ENST00000287936 usage") + scale_y_continuous(limits = c(0.4,0.9))
ggsave("results/figures/HMGCR_ENST00000287936_ensembl.pdf", plot = boxplot, width = 2.5, height = 3)




#Extract HMGCR phenotype
selected_phenotype_id = "ENSG00000113161.grp_2.contained.ENST00000343975"
selected_snp_id = "rs3846662"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N"))

boxplot = plotQtlCol(plot_data) + ylab("txrevse: ENST00000343975 usage") + scale_y_continuous(limits = c(0,0.5))
ggsave("results/figures/HMGCR_ENST00000343975_txrevise.pdf", plot = boxplot, width = 2.5, height = 3)



#Make a boxplot of the other trancript
selected_phenotype_id2 = "ENSG00000113161.grp_2.contained.ENST00000511206"

plot_data = constructQtlPlotDataFrame(selected_phenotype_id2, selected_snp_id, 
                                      ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N"))

boxplot2 = plotQtlCol(plot_data) + ylab("txrevise: ENST00000511206 usage")
ggsave("results/figures/HMGCR_boxplot2.pdf", plot = boxplot2, width = 2.5, height = 3)


#Plot transcripts
hmgcr_transcripts = dplyr::filter(ensembl_gene_meta, gene_name == "HMGCR")
hmgcr_tx_plot = plotTranscripts(exons[hmgcr_transcripts$transcript_id],cdss[intersect(names(cdss), hmgcr_transcripts$transcript_id)])
ggsave("results/figures/HMGCR_transcripts.pdf", plot = hmgcr_tx_plot, width = 6, height = 5)


