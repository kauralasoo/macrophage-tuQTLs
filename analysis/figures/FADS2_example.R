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
se_featureCounts = readRDS("results/SummarizedExperiments/acLDL_featureCounts.rds")
cqn_matrix = assays(se_featureCounts)$cqn
gene_meta = rowData(se_featureCounts) %>% tbl_df2()
sample_meta = colData(se_featureCounts) %>% tbl_df2()

#Import genotypes
vcf_file = readRDS("results/genotypes/acLDL/imputed.70_samples.sorted.filtered.named.rds")
GRCh38_information = importVariantInformation("results/genotypes/acLDL/imputed.70_samples.variant_information.txt.gz")
GRCh37_information = importVariantInformation("results/genotypes/acLDL/GRCh37/imputed.70_samples.variant_information.GRCh37.txt.gz")

#Import coloc resuts
coloc_hits = readRDS("results/coloc/acLDL_GWAS_coloc_hits.rds")

#Extract FADS2
phenotype_id = "ENSG00000134824"
snp_id = "rs968567"

#Make a QTL boxplot
plot_data = constructQtlPlotDataFrame(phenotype_id, snp_id, 
                                      cqn_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name)

boxplot = plotQtlCol(plot_data)
ggsave("results/figures/FADS2_acLDL_boxplot.pdf", plot = boxplot, width = 2.5, height = 4)



#### Make manhattan plots for the QTL
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name", "type"))
region_coords = c(61629296, 62027885)

#Import QTL and GWAS summary stats and convert them to the same GRCh38 coordinate space

##### Promoter QTL ####
qtl_df = data_frame(phenotype_id = phenotype_id, snp_id = snp_id, trait = "RA")
qtl_paths = list(Ctrl = "processed/acLDL/qtltools/output/featureCounts/sorted/Ctrl.nominal.sorted.txt.gz",
                 AcLDL = "processed/acLDL/qtltools/output/featureCounts/sorted/AcLDL.nominal.sorted.txt.gz")
qtl_summary = importSummariesForPlotting(qtl_df, gwas_stats_labeled, 
                                         gwas_dir = "/Users/alasoo/datasets/Inflammatory_GWAS/", qtl_paths = qtl_paths, 
                                         GRCh37_variants = GRCh37_information, GRCh38_variants = GRCh38_information, 
                                         cis_dist = 2e5, type = "QTLTools") %>%
  arrange(condition_name, p_nominal) %>% 
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::mutate(track_id = condition_name)

#Make a manhattan plot
gwas_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id == "RA"), region_coords, color_R2 = TRUE)
eqtl_manhattan = wiggleplotr::makeManhattanPlot(dplyr::filter(qtl_summary, track_id != "RA"), region_coords, color_R2 = TRUE, data_track = FALSE) +
  theme(plot.margin=unit(c(0.1,1,0.1,1),"line"),
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(colour = "grey10"),
      strip.background = element_rect(fill = "grey85")) +
  xlab("Chromosome 11 position")


joint_plot = cowplot::plot_grid(gwas_manhattan, eqtl_manhattan, 
                                align = "v", ncol = 1, rel_heights = c(1,2))
ggsave("results/figures/FADS2_acLDL_manhattan.pdf", plot = joint_plot, width = 4, height = 6)


