library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("data.table")
library("GenomicFeatures")
load_all("../seqUtils/")
load_all("../wiggleplotr")

#Import SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")

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




