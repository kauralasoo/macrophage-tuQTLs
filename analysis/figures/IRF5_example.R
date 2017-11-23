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
GRCh38_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
GRCh37_information = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")

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
tx_metadata = gene_meta %>% tbl_df2() %>% 
  dplyr::filter(gene_name == "IRF5") %>% 
  tidyr::separate(group_id, c("gene","group","position"), "\\.") %>%
  dplyr::select(-gene) %>%
  dplyr::filter(group == "grp_1")

#Upstream
txs = dplyr::filter(tx_metadata, position == "upstream")
tx_ranges = revised_granges[txs$transcript_id]
tx_plot = plotTranscripts(tx_ranges, transcript_label = FALSE)
ggsave("results/figures/IRF5_promoter_tx.pdf", plot = tx_plot, width = 3, height = 1.5)

#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs3778754", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(track_id = figure_name)

#Make a zoomed coverage plot
redundant_exons = purrr::reduce(tx_ranges, intersect)[1:7]
new_ranges = purrr::map(as.list(tx_ranges), ~setdiff(.,redundant_exons))
plotTranscripts(new_ranges[c(2,5)], rescale_introns = FALSE)

coverage_plot = plotCoverage(exons = new_ranges[c(2,5)], 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.5,0.5), 
                             return_subplots_list = TRUE, coverage_type = "line", transcript_label = FALSE)
ggsave("results/figures/IRF5_promoter_coverage.pdf", plot = coverage_plot$coverage_plot, width = 3, height = 1.2)


#Make a boxplot
boxplot = constructQtlPlotDataFrame("ENSG00000128604.grp_1.upstream.ENST00000249375", "rs3778754", 
                                            tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(constructGenotypeText("rs3778754", variant_information), by = "genotype_value") %>%
  plotQtlCol() + 
  ylab("Promoter splicing")
ggsave("results/figures/IRF5_promoter_boxplot.pdf", plot = boxplot, width = 2.5, height = 2.5)


###### Contained ######
txs = dplyr::filter(tx_metadata, position == "contained")
tx_ranges = revised_granges[txs$transcript_id]
tx_plot = plotTranscripts(tx_ranges, transcript_label = FALSE)
ggsave("results/figures/IRF5_RI_tx.pdf", plot = tx_plot, width = 3, height = 1.5)

#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs199508964", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(track_id = figure_name)

#Select relevant exons
new_ranges = list(tx1 = tx_ranges[[7]][5], tx2 = tx_ranges[[8]][5:6])

#Make coverage plot
coverage_plot = plotCoverage(exons = new_ranges, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.5,0.5), 
                             return_subplots_list = TRUE, coverage_type = "line", transcript_label = FALSE)
ggsave("results/figures/IRF5_RI_coverage.pdf", plot = coverage_plot$coverage_plot, width = 3, height = 1.2)

#Make a boxplot
boxplot = constructQtlPlotDataFrame("ENSG00000128604.grp_1.contained.ENST00000619830", "rs199508964", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(constructGenotypeText("rs199508964", variant_information), by = "genotype_value") %>%
  dplyr::mutate(genotype_text = genotype_value) %>%
  plotQtlCol() + 
  ylab("Intron retention")
ggsave("results/figures/IRF5_RI_boxplot.pdf", plot = boxplot, width = 2.5, height = 2.5)


###### Downstream ######
txs = dplyr::filter(tx_metadata, position == "downstream")
tx_ranges = revised_granges[txs$transcript_id]
tx_plot = plotTranscripts(tx_ranges, transcript_label = FALSE)
ggsave("results/figures/IRF5_UTR_tx.pdf", plot = tx_plot, width = 3, height = 1.5)

#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str2_df, "rs10954213", vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive"))

redundant_exons = purrr::reduce(tx_ranges, intersect)
new_ranges = purrr::map(as.list(tx_ranges), ~setdiff(.,redundant_exons))[2:4]
plotTranscripts(new_ranges)

#Make coverage plot
coverage_plot = plotCoverage(exons = new_ranges, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.5,0.5), 
                             return_subplots_list = TRUE, coverage_type = "line", transcript_label = FALSE)
ggsave("results/figures/IRF5_UTR_coverage.pdf", plot = coverage_plot$coverage_plot, width = 3, height = 1.2)

#Make a boxplot
boxplot = constructQtlPlotDataFrame("ENSG00000128604.grp_1.downstream.ENST00000489702", "rs10954213", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::filter(condition_name %in% c("naive")) %>%
  dplyr::left_join(constructGenotypeText("rs10954213", variant_information), by = "genotype_value") %>%
  plotQtlCol() +
  ylab("Alternative 3'UTR")
ggsave("results/figures/IRF5_UTR_boxplot.pdf", plot = boxplot, width = 2.5, height = 2.5)



#### Make manhattan plots for each QTL
#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name", "type"))
region_coords = c(128749500, 129149286)

#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space

##### Promoter QTL ####
qtl_df = data_frame(phenotype_id = "ENSG00000128604.grp_1.upstream.ENST00000249375", snp_id = "rs3778754", trait = "RA")
qtl_paths = list(naive = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/reviseAnnotations/sorted/naive.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, 
                                         gwas_dir = "~/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = GRCh37_information, GRCh38_variants = GRCh38_information, 
                                         cis_dist = 2e5, type = "QTLTools") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = as.character(track_id)) %>%
  dplyr::mutate(track_id = ifelse(track_id == "naive", "promoter splicing", track_id)) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("RA","promoter splicing")))

#Make a manhattan plot
gwas_lead_var = dplyr::filter(qtl_summary, snp_id == "rs3778753", track_id == "RA")
gwas_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id == "RA"), region_coords, color_R2 = TRUE) +
  geom_point(data = gwas_lead_var, color = "red")

qtl_lead_var = dplyr::filter(qtl_summary, snp_id == "rs3778753", track_id != "RA")
promoter_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id != "RA"), region_coords, color_R2 = TRUE) + 
  geom_point(data = qtl_lead_var, color = "red")


##### Retained intron QTL #####
qtl_df = data_frame(phenotype_id = "ENSG00000128604.grp_1.contained.ENST00000619830", snp_id = "rs199508964", trait = "RA")
qtl_paths = list(naive = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/reviseAnnotations/sorted/naive.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, 
                                         gwas_dir = "~/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = GRCh37_information, GRCh38_variants = GRCh38_information, 
                                         cis_dist = 2e5, type = "QTLTools") %>%
  dplyr::mutate(condition_name = as.character(condition_name)) %>%
  dplyr::arrange(desc(condition_name), p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = ifelse(track_id == "naive", "retained intron", track_id)) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("RA","retained intron")))

#Make a manhattan plot
qtl_lead_var = dplyr::filter(qtl_summary, snp_id == "rs3778753", track_id != "RA")
RI_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id != "RA"), region_coords, color_R2 = TRUE) + 
  geom_point(data = qtl_lead_var, color = "red")



##### 3' UTR #####
qtl_df = data_frame(phenotype_id = "ENSG00000128604.grp_1.downstream.ENST00000489702", snp_id = "rs10954213", trait = "RA")
qtl_paths = list(naive = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/reviseAnnotations/sorted/naive.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, 
                                         gwas_dir = "~/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = GRCh37_information, GRCh38_variants = GRCh38_information, 
                                         cis_dist = 2e5, type = "QTLTools") %>%
  dplyr::mutate(condition_name = as.character(condition_name)) %>%
  dplyr::arrange(desc(condition_name), p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::rename(track_id = condition_name) %>%
  dplyr::mutate(track_id = ifelse(track_id == "naive", "alt 3'UTR", track_id)) %>%
  dplyr::mutate(track_id = factor(track_id, levels = c("RA","alt 3'UTR")))

#Make a manhattan plot
qtl_lead_var = dplyr::filter(qtl_summary, snp_id == "rs3778753", track_id != "RA")
UTR_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id != "RA"), region_coords, color_R2 = TRUE, data_track = FALSE) + 
  geom_point(data = qtl_lead_var, color = "red") +
  theme(plot.margin=unit(c(0.1,1,0.1,1),"line"),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(colour = "grey10"),
        strip.background = element_rect(fill = "grey85")) +
  xlab("Chromosome 7 position")


joint_plot = cowplot::plot_grid(gwas_manhattan, promoter_manhattan, RI_manhattan, UTR_manhattan, 
                                align = "v", ncol = 1, rel_heights = c(1,1,1,1.5))
ggsave("results/figures/IRF5_manhattan_plot.pdf", width = 4, height = 8)



#Import coloc posterior probabilities:
coloc_hits = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")
View(dplyr::filter(coloc_hits$reviseAnnotations, gene_name == "IRF5"))
