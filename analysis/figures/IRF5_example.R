library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("../wiggleplotr")

#Import SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")

#Extract TPM matrix
tpm_matrix = assays(se_revised)$tpm_ratios
sample_meta = colData(se_revised) %>% tbl_df2()
gene_meta = rowData(se_revised) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = transcript_id)

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")

#Set up bigwig metadata
rna_meta_str1_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                                tbl_df2(colData(se_featureCounts)),
                                                "/Volumes/Ajamasin/bigwig/salmonella", 
                                                bigWig_suffix = ".str1.bw")
rna_meta_str2_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                                tbl_df2(colData(se_featureCounts)),
                                                "/Volumes/Ajamasin/bigwig/salmonella", 
                                                bigWig_suffix = ".str2.bw")

#Import GRangesList
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")

#Extract transcript metadata
tx_metadata = rowData(se_revised) %>% tbl_df2() %>% 
  dplyr::filter(gene_name == "IRF5") %>% 
  tidyr::separate(group_id, c("gene","group","position"), "\\.") %>%
  dplyr::select(-gene) %>%
  dplyr::filter(group == "grp_1")

#Upstream
txs = dplyr::filter(tx_metadata, position == "upstream")
tx_ranges = revised_granges[txs$transcript_id]
plotTranscripts(tx_ranges)

#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs3778754", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive"))

#Make coverage plot
coverage_plot = plotCoverage(exons = tx_ranges, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.55,0.45), 
                             return_subplots_list = FALSE, coverage_type = "line")
ggsave("results/figures/IRF5_promoter_coverage.pdf", plot = coverage_plot, width = 6, height = 7)

#Make a boxplot
boxplot = constructQtlPlotDataFrame("ENSG00000128604.grp_1.upstream.ENST00000249375", "rs3778754", 
                                            tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(constructGenotypeText("rs3778754", variant_information), by = "genotype_value") %>%
  plotQtlCol()
ggsave("results/figures/IRF5_promoter_boxplot.pdf", plot = boxplot, width = 4, height = 4)


###### Contained ######
txs = dplyr::filter(tx_metadata, position == "contained")
tx_ranges = revised_granges[txs$transcript_id]
plotTranscripts(tx_ranges)

#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs199508964", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive"))

#Make coverage plot
coverage_plot = plotCoverage(exons = tx_ranges, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.55,0.45), 
                             return_subplots_list = FALSE, coverage_type = "line")
ggsave("results/figures/IRF5_RI_coverage.pdf", plot = coverage_plot, width = 6, height = 7)

#Make a boxplot
boxplot = constructQtlPlotDataFrame("ENSG00000128604.grp_1.contained.ENST00000619830", "rs199508964", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(constructGenotypeText("rs199508964", variant_information), by = "genotype_value") %>%
  plotQtlCol()
ggsave("results/figures/IRF5_RI_boxplot.pdf", plot = boxplot, width = 4, height = 4)


###### Downstream ######
txs = dplyr::filter(tx_metadata, position == "downstream")
tx_ranges = revised_granges[txs$transcript_id]
plotTranscripts(tx_ranges)

#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs10954213", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive"))

#Make coverage plot
coverage_plot = plotCoverage(exons = tx_ranges, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.55,0.45), 
                             return_subplots_list = FALSE, coverage_type = "line")
ggsave("results/figures/IRF5_UTR_coverage.pdf", plot = coverage_plot, width = 6, height = 5)

#Make a boxplot
boxplot = constructQtlPlotDataFrame("ENSG00000128604.grp_1.downstream.ENST00000489702", "rs10954213", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(constructGenotypeText("rs10954213", variant_information), by = "genotype_value") %>%
  plotQtlCol()
ggsave("results/figures/IRF5_UTR_boxplot.pdf", plot = boxplot, width = 4, height = 4)



