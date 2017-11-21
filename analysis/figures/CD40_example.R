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

#Import transcript coordinates
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names=TRUE)
cdss = cdsBy(txdb, by = "tx", use.names=TRUE)

#Extract TPM matrix
ratio_matrix = assays(se_revised)$tpm_ratios
tpm_matrix = assays(se_revised)$tpms
sample_meta = colData(se_revised) %>% tbl_df2()
gene_meta = rowData(se_revised) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = transcript_id)

#Import coloc resuts
coloc_hits = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
GRCh38_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
GRCh37_information = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")

#Extract CD40
selected_phenotype_id = "ENSG00000101017.grp_1.upstream.ENST00000372285"
selected_snp_id = "rs4239702"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N","I"))

boxplot = plotQtlCol(plot_data) +
  ylab("txrevise: 372285")
ggsave("results/figures/CD40_promoter_relative_expression.pdf", boxplot, width = 1.75, height = 2.5)


#Make a QTL boxplot (absolute expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N","I"))

boxplot = plotQtlCol(plot_data) +
  ylab("ENST00000372285 promoter expression (TPM)")
ggsave("results/figures/CD40_ENST00000372285_absolute_expression.pdf", boxplot, width = 2.5, height = 3.5)

#The other promoter
plot_data = constructQtlPlotDataFrame("ENSG00000101017.grp_1.upstream.ENST00000372276", selected_snp_id, 
                                      tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N","I"))

boxplot = plotQtlCol(plot_data) +
  ylab("ENST00000372276 promoter expression (TPM)")
ggsave("results/figures/CD40_ENST00000372276_absolute_expression.pdf", boxplot, width = 2.5, height = 3.5)



#Make a QTL boxplot (eQTL)
plot_data = constructQtlPlotDataFrame("ENSG00000101017", selected_snp_id, 
                                      assays(se_featureCounts)$cqn, vcf_file$genotypes, sample_meta, rowData(se_featureCounts) %>% tbl_df2()) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name) %>%
  dplyr::filter(condition_name %in% c("N","I"))

boxplot = plotQtlCol(plot_data)
ggsave("results/figures/CD40_eQTL_boxplot.pdf", boxplot, width = 2.5, height = 3.5)



#Make a manhattan plot
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name", "type"))

##### Promoter QTL ####
qtl_df = data_frame(phenotype_id = selected_phenotype_id, snp_id = selected_snp_id, trait = "RA")
qtl_paths = list(naive = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/reviseAnnotations/sorted/naive.nominal.sorted.txt.gz",
                 IFNg = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/reviseAnnotations/sorted/IFNg.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, 
                                         gwas_dir = "/Users/alasoo/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = GRCh37_information, GRCh38_variants = GRCh38_information, 
                                         cis_dist = 2e5, type = "QTLTools") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>% 
  dplyr::mutate(track_id = condition_name)

#Make a manhattan plot
region_coords = c(45980000, 46200000)
gwas_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id == "RA"), region_coords, color_R2 = TRUE)
eqtl_data = dplyr::filter(qtl_summary, track_id != "RA") %>% 
  dplyr::left_join(conditionFriendlyNames()) %>% 
  dplyr::mutate(track_id = figure_name)
eqtl_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(eqtl_data, track_id != "RA"), region_coords, color_R2 = TRUE)

#Make a gene location plot
cd40_exons = ensembl_exons[1]
names(cd40_exons) = c("CD40")
transcripts_plot = plotTranscripts(exons = cd40_exons, rescale_introns = FALSE, region_coords = region_coords)

joint_plot = cowplot::plot_grid(gwas_manhattan, eqtl_manhattan, transcripts_plot,
                                align = "v", ncol = 1, rel_heights = c(3,6,3))


ggsave("results/figures/CD40_manhattan.pdf", plot = joint_plot, width = 4, height = 5)



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
  dplyr::filter(gene_name == "CD40") %>%
  dplyr::filter(transcript_biotype == "protein_coding") %>%
  dplyr::filter(transcript_id %in% c("ENST00000372285","ENST00000372276"))
ensembl_exons = exons[ensembl_meta$transcript_id] %>% txrevise::removeMetadata()
ensembl_cdss = cdss[ensembl_meta$transcript_id] %>% txrevise::removeMetadata()
names(ensembl_exons) = paste("CD40:", names(ensembl_exons), sep = "")
names(ensembl_cdss) = paste("CD40:", names(ensembl_cdss), sep = "")

##### txrevise #####
#Extract transcript metadata
tx_metadata = gene_meta %>% tbl_df2() %>% 
  dplyr::filter(gene_name == "CD40") %>% 
  tidyr::separate(group_id, c("gene","group","position"), "\\.") %>%
  dplyr::select(-gene) %>%
  dplyr::filter(group == "grp_1") %>%
  dplyr::filter(position == "upstream")

tx_ranges = revised_granges[tx_metadata$transcript_id]
names(tx_ranges) = c("txrevise: ENST00000372285", "txrevise: ENST00000372276")


#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, selected_snp_id, vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive", "IFNg")) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(track_id = figure_name)

#Make coverage plot
transcripts_plot = plotTranscripts(exons = c(ensembl_exons, tx_ranges), cdss = ensembl_cdss)
ggsave("results/figures/CD40_transcripts.pdf", transcripts_plot, width = 3.5, height = 2.5)

#Make a zoomed-in coverage plot of the first two exons
first_exons = list(tx_ranges[[1]][1:2], tx_ranges[[2]][1:2])
names(first_exons) = names(tx_ranges)

coverage_plot = plotCoverage(exons = first_exons, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.7,0.3), 
                             return_subplots_list = F, coverage_type = "line", transcript_label = FALSE)
ggsave("results/figures/CD40_coverage_zoomed.pdf", coverage_plot, width = 3, height = 3)



