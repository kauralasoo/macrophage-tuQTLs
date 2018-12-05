library("dplyr")
library("readr")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../txrevise/")
load_all("../seqUtils/")
load_all("../wiggleplotr/")
load_all("~/software/rasqual/rasqualTools/")

#Import overlapping promoters
colocalised_pairs = read_tsv("results/tables/colocalised_puQTL_caQTL_pairs.txt")
qtl_hits = dplyr::transmute(colocalised_pairs, phenotype_id = transcript_id, snp_id, group_id)

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")

#Import Ensembl transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Import revisedAnnotations Granges
revised_granges = updateObject(readRDS("results/annotations/txrevise_promoters.GRangesList.rds"), verbose = T)
leafcutter_granges = readRDS("results/annotations/salmonella_leafcutter.GRangesList.rds")

#Import gene exrression data
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
revised_gene_meta = rowData(se_revised) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = group_id)
revised_groups = dplyr::transmute(revised_gene_meta, phenotype_id = transcript_id, group_id)


#Construct dfs with paths to bigwig files
str1_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                       colData(se_revised) %>% tbl_df2(), 
                                       "/Volumes/Ajamasin/bigwig/salmonella/", 
                                       bigWig_suffix = ".str1.bw",
                                       condition_name_levels = c("naive","IFNg", "SL1344", "IFNg_SL1344")) %>%
  dplyr::filter(condition_name %in% c("naive"))
str2_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                       colData(se_revised) %>% tbl_df2(), 
                                       "/Volumes/Ajamasin/bigwig/salmonella/", 
                                       bigWig_suffix = ".str2.bw",
                                       condition_name_levels = c("naive","IFNg", "SL1344", "IFNg_SL1344")) %>%
  dplyr::filter(condition_name %in% c("naive"))

#Make plots
selected_hit = dplyr::filter(qtl_hits, phenotype_id == "ENSG00000115677.grp_2.upstream.ENST00000391976")[1,]
plots_df = purrrlyr::by_row(selected_hit, ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes,
                                                          revised_granges, 
                                                          gene_metadata = revised_gene_meta,
                                                          plot_fraction = 0.2, coverage_type = "line", 
                                                          rescale_introns = TRUE, heights = c(0.5,0.5)))
plots_list = plots_df$.out
names(plots_list) = plots_df$phenotype_id
savePlotList(plots_list, "processed/salmonella/coverage_plots/puQTL_caQTL_olaps/", suffix = ".pdf", width = 10, height = 10)



#Visualise caQTLs from the same loci
#Import ATAC data
atac_list = readRDS("../macrophage-gxe-study/results/ATAC/ATAC_combined_accessibility_data.rds")
atac_list$sample_metadata$condition_name = factor(atac_list$sample_metadata$condition_name, 
                                                  levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344"))
atac_meta_df = wiggleplotrConstructMetadata(atac_list$counts, atac_list$sample_metadata, "/Volumes/Ajamasin/bigwig/ATAC/")

#Define region of interest
region_coords = c(241314000,241317000)

#Fetch all peaks in the region and make peak annot plot
peak_annot = wiggleplotrExtractPeaks(region_coords, chrom = 2, atac_list$gene_metadata)
peak_plot = plotTranscripts(peak_annot$peak_list, peak_annot$peak_list, peak_annot$peak_annot, rescale_introns = FALSE, 
                            region_coords = region_coords, connect_exons = FALSE, transcript_label = FALSE) + dataTrackTheme()

#Construct metadata df for wiggleplotr
atac_track_data = wiggleplotrGenotypeColourGroup(atac_meta_df, "rs12624195", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive")) %>%
  dplyr::mutate(track_id = "ATAC-seq")

ATAC_coverage = plotCoverage(exons = peak_annot$peak_list, cdss = peak_annot$peak_list, track_data = atac_track_data, rescale_introns = FALSE, 
                             transcript_annotations = peak_annot$peak_annot, fill_palette = getGenotypePalette(),
                             connect_exons = FALSE, transcript_label = FALSE, plot_fraction = 0.1, heights = c(0.7,0.3), 
                             region_coords = region_coords, return_subplots_list = TRUE, coverage_type = "line")


#HDLBP
#Plot all tx of this gene
#Make plots
selected_hit = dplyr::filter(qtl_hits, phenotype_id == "ENSG00000115677.grp_2.upstream.ENST00000391976")[1,]
dplyr::mutate(selected_hit, snp_id =  "rs12624195")
plots_df = purrrlyr::by_row(selected_hit, ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes,
                                                               revised_granges, transcript_label = FALSE,
                                                               gene_metadata = revised_gene_meta,
                                                               plot_fraction = 0.2, coverage_type = "line", 
                                                               rescale_introns = TRUE, heights = c(0.4,0.6)))
plots_list = plots_df$.out
names(plots_list) = plots_df$phenotype_id
ggsave("results/figures/HDLBP_coverage.pdf", plot = plots_list[[1]], width = 5, height = 4)

tx_meta = dplyr::filter(revised_gene_meta, group_id == "ENSG00000115677.grp_2.upstream")
transcripts = as.list(revised_granges[tx_meta$transcript_id])
shared_exons = listIntersect(transcripts)
shared_exons = sort(shared_exons)
promoters = purrr::map(transcripts, ~setdiff(., shared_exons))

#keep only relevant
selected_ids = c("ENSG00000115677.grp_2.upstream.ENST00000441124", "ENSG00000115677.grp_2.upstream.ENST00000391976",
                 "ENSG00000115677.grp_2.upstream.ENST00000428482","ENSG00000115677.grp_2.upstream.ENST00000420451",
                 "ENSG00000115677.grp_2.upstream.ENST00000310931")
promoters = promoters[selected_ids]

promoter_ATAC_coverage = plotCoverage(exons = promoters, cdss = promoters, track_data = atac_track_data, rescale_introns = FALSE, 
                             fill_palette = getGenotypePalette(),
                             connect_exons = TRUE, transcript_label = FALSE, plot_fraction = 0.1, heights = c(0.7,0.3), 
                             region_coords = region_coords,
                             return_subplots_list = TRUE, coverage_type = "line")


#Make promoter RNA coverage plot
rna_track_data = wiggleplotrGenotypeColourGroup(str1_df, "rs12624195", vcf_file$genotypes, 1) %>%
  dplyr::filter(track_id %in% c("naive")) %>%
  dplyr::mutate(track_id = "RNA-seq")
promoter_RNA_plot = plotCoverage(exons = promoters, cdss = promoters, track_data = rna_track_data, rescale_introns = FALSE, 
                                 fill_palette = getGenotypePalette(),
                                 connect_exons = TRUE, transcript_label = FALSE, plot_fraction = 0.1, heights = c(0.7,0.3), 
                                 region_coords = region_coords,
                                 return_subplots_list = TRUE, coverage_type = "line")


#Make manhattan plot of the caQTL p-values
atac_summary = tabixFetchGenesQuick(c("ATAC_peak_156061"),
                      tabix_file = "~/databases/ATAC/rasqual/naive_100kb.sorted.txt.gz", 
                      gene_metadata = atac_list$gene_metadata, cis_window = 1e5)[[1]] %>%
  dplyr::mutate(track_id = "ATAC-seq") %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes) 
region_coords2 = c(min(atac_summary$pos), max(atac_summary$pos))
peak_manhattan = makeManhattanPlot(atac_summary, region_coords, color_R2 = TRUE)

#Make a joint plot
joint_plot = cowplot::plot_grid(peak_manhattan, ATAC_coverage$coverage_plot, promoter_RNA_plot$coverage_plot, ATAC_coverage$tx_structure + dataTrackTheme(), promoter_ATAC_coverage$tx_structure, 
                                align = "v", ncol = 1, rel_heights = c(3,3,3,1,3))
ggsave("results/figures/HDLBP_promoter_accessibility_plot.pdf", plot = joint_plot, width = 3.5, height = 6)


#Make a manhattan plot of the puQTL
selected_hit = dplyr::filter(qtl_hits, phenotype_id == "ENSG00000115677.grp_2.upstream.ENST00000391976")[1,] 
qtl_ranges = constructVariantRanges(selected_hit, variant_information, cis_dist = 1e6)
qtl_path = "/Volumes/Ajamasin/processed/salmonella/qtltools/output/txrevise_promoters/sorted/naive.nominal.sorted.txt.gz"
qtl_pvalues = qtltoolsTabixFetchPhenotypes(qtl_ranges, qtl_path)[[1]] %>%
  dplyr::mutate(track_id = "HDLBP puQTL") %>%
  dplyr::arrange(p_nominal) %>%
  addR2FromLead(vcf_file$genotypes) %>%
  dplyr::mutate(pos = snp_start) %>%
  dplyr::mutate(lead_variant = ifelse(snp_id == "rs62187434", "puQTL", NA)) %>%
  dplyr::mutate(lead_variant = ifelse(snp_id == "rs12624195", "caQTL", lead_variant))
region_coords2 = c(min(qtl_pvalues$pos), max(qtl_pvalues$pos))
region_coords2 = c(241000000,241700000)
peak_manhattan = makeManhattanPlot(qtl_pvalues, region_coords2, color_R2 = TRUE, data_track = FALSE) +
  theme(strip.text.y = element_text(colour = "grey10"),
        strip.background = element_rect(fill = "grey85"))

#Add lead vars
lead_vars = dplyr::filter(qtl_pvalues, !is.na(lead_variant))
updated_manhattan = peak_manhattan + geom_point(data = lead_vars, color = "red")
ggsave("results/figures/HDLBP_manhattan.pdf", plot = updated_manhattan, width = 5, height = 2)


calculatePairR2("rs62187434", "rs12624195", vcf_file$genotypes)


#Make caQTL boxplot
plot_data = constructQtlPlotDataFrame("ATAC_peak_156061", "rs12624195", 
                                      atac_list$cqn, vcf_file$genotypes, atac_list$sample_metadata, atac_list$gene_metadata) %>% 
  dplyr::left_join(constructGenotypeText("rs12624195", variant_information), by = "genotype_value") %>%
  dplyr::filter(condition_name == "naive")
plotQtlCol(plot_data)

