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

fc_gene_meta = rowData(se_featureCounts) %>% tbl_df2()
fc_sample_meta = colData(se_featureCounts) %>% tbl_df2()
cqn_matrix = assays(se_featureCounts)$cqn

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

#Visualise some of them
#Extract HAUS8 phenotype
selected_phenotype_id = "ENSG00000131351.grp_1.upstream.ENST00000597917"
selected_snp_id = "rs146734736"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      ratio_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name)  %>%
  dplyr::filter(condition_name %in% c("N", "I+S"))

boxplot = plotQtlCol(plot_data) + scale_color_manual(values = conditionPalette()[c(1,4)]) +
  theme(legend.position = "None") +
  ylab("ENST00000597917 promoter usage")
ggsave("results/figures/HAUS8_tpm_ratio.pdf", plot = boxplot, width = 2, height = 3)


#Make a QTL boxplot (absolute expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name)  %>%
  dplyr::filter(condition_name %in% c("N", "I+S"))

boxplot = plotQtlCol(plot_data) + scale_color_manual(values = conditionPalette()[c(1,4)]) +
  theme(legend.position = "None") +
  ylab("ENST00000597917 promoter expression")
ggsave("results/figures/HAUS8_tpm.pdf", plot = boxplot, width = 2, height = 3)


#Make a read coverage plot
txrevise_ids = dplyr::filter(gene_meta, group_id == "ENSG00000131351.grp_1.upstream")$transcript_id
ensembl_ids = filter(gene_meta_ensembl, gene_id == "ENSG00000131351")$transcript_id

selected_txrevise = c("ENSG00000131351.grp_1.upstream.ENST00000597917", "ENSG00000131351.grp_1.upstream.ENST00000253669")
plotTranscripts(revised_granges[selected_txrevise])
plotTranscripts(exons[ensembl_ids])


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
  dplyr::filter(gene_name == "HAUS8") %>%
  dplyr::filter(transcript_biotype == "protein_coding") %>%
  dplyr::filter(transcript_id %in% c("ENST00000253669", "ENST00000593360"))
ensembl_exons = exons[ensembl_meta$transcript_id] %>% txrevise::removeMetadata()
ensembl_cdss = cdss[ensembl_meta$transcript_id] %>% txrevise::removeMetadata()
names(ensembl_exons) = paste("HAUS8: ", names(ensembl_exons), sep = "")
names(ensembl_cdss) = paste("HAUS8: ", names(ensembl_cdss), sep = "")

tx_ranges = revised_granges[selected_txrevise]
names(tx_ranges) = c("txrevise: ENST00000597917", "txrevise: ENST00000253669")


#Add genotype data
track_data = wiggleplotrGenotypeColourGroup(rna_meta_str1_df, selected_snp_id, vcf_file$genotypes, -1)
filtered_tracks = dplyr::filter(track_data) %>% dplyr::filter(condition_name %in% c("naive", "IFNg_SL1344")) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(track_id = figure_name)


coverage_plot = plotCoverage(exons = c(ensembl_exons, rev(tx_ranges)), cdss = ensembl_cdss, 
                             track_data = filtered_tracks, rescale_introns = TRUE,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.2, heights = c(0.5,0.5), 
                             return_subplots_list = F, coverage_type = "line", transcript_label = TRUE)
ggsave("results/figures/HAUS8_coverage_plot.pdf", plot = coverage_plot, width = 6, height = 4)



#MAke eQTL boxplot
selected_phenotype_id = "ENSG00000131351"
selected_snp_id = "rs146734736"

#Make a QTL boxplot (relative expression)
plot_data = constructQtlPlotDataFrame(selected_phenotype_id, selected_snp_id, 
                                      cqn_matrix, vcf_file$genotypes, fc_sample_meta, fc_gene_meta) %>% 
  dplyr::left_join(constructGenotypeText(selected_snp_id, GRCh38_information), by = "genotype_value") %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::mutate(condition_name = figure_name)  %>%
  dplyr::filter(condition_name %in% c("N", "I+S"))

boxplot = plotQtlCol(plot_data) + scale_color_manual(values = conditionPalette()[c(1,4)]) +
  theme(legend.position = "None") +
  ylab("HAUS8 expression")
ggsave("results/figures/HAUS8_gene_expression.pdf", plot = boxplot, width = 2, height = 3)

