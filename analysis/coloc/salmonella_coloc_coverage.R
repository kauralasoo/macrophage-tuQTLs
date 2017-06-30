library("dplyr")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

#Import gene expression data
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")

#Import SummarizedExperiments for all phenotypes
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
se_reviseAnnotations = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
se_list = list(Ensembl_87 = se_ensembl, reviseAnnotation = se_reviseAnnotations, leafcutter = se_leafcutter)

#Import QTL mapping results
trqtl_min_pvalues = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")[1:3]

#Extract sample metadata
sample_meta_list = purrr::map(se_list, ~colData(.) %>% tbl_df2())
gene_meta_list = purrr::map(se_list, ~rowData(.) %>% tbl_df2())
gene_meta_list$reviseAnnotation = dplyr::mutate(gene_meta_list$reviseAnnotation, gene_id = group_id)

#Import Ensembl transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Import revised and leafcutter annotations
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")
leafcutter_granges = readRDS("results/annotations/salmonella_leafcutter.GRangesList.rds")

#Set up bigwig metadata
str1_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                                tbl_df2(colData(se_featureCounts)),
                                                "/Volumes/Ajamasin/bigwig/salmonella", 
                                                bigWig_suffix = ".str1.bw")
str2_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                                tbl_df2(colData(se_featureCounts)),
                                                "/Volumes/Ajamasin/bigwig/salmonella", 
                                                bigWig_suffix = ".str2.bw")

#Import GWAS overlaps
gwas_olaps = readRDS("results/coloc/salmonella_GWAS_coloc_hits.rds")

#Visualise leafCutter junctions
leafcutter_olaps = dplyr::select(gwas_olaps$leafcutter, phenotype_id, snp_id, gene_name) %>% unique()
leafcutter_hits = dplyr::semi_join(purrr::map_df(trqtl_min_pvalues$leafcutter, identity), leafcutter_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::left_join(leafcutter_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::select(group_id,phenotype_id, snp_id, gene_name) %>%
  unique()

#Make all coverage plots
plots_df = purrrlyr::by_row(leafcutter_hits[1:2,], ~makeQTLCoveragePlot(.,str1_df, str2_df, vcf_file$genotypes,
                                                               leafcutter_granges, 
                                                               gene_metadata = gene_meta_list$leafcutter,
                                                               plot_fraction = 0.2, coverage_type = "line", 
                                                               rescale_introns = TRUE, heights = c(0.6,0.4)), .to = "plot")
plots_df_df = dplyr::mutate(plots_df, plot_title = paste(gene_name, snp_id, phenotype_id, sep = "_"))
plot_list = setNames(plots_df_df$plot, plots_df_df$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/coverage/leafcutter/")

#Visualise reviseAnnotations
revised_olaps = dplyr::select(gwas_olaps$revisedAnnotation, phenotype_id, snp_id, gene_name) %>% unique()
revised_hits = dplyr::semi_join(purrr::map_df(trqtl_min_pvalues$reviseAnnotations, identity), revised_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::left_join(revised_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::select(phenotype_id, snp_id, gene_name) %>%
  dplyr::left_join(dplyr::transmute(gene_meta_list$reviseAnnotation, phenotype_id = transcript_id, group_id), by = "phenotype_id") %>%
  unique()

#Make all coverage plots
plots_df = purrrlyr::by_row(revised_hits, ~makeQTLCoveragePlot(.,str1_df, str2_df, vcf_file$genotypes,
                                                            revised_granges, 
                                                            gene_metadata = gene_meta_list$reviseAnnotation,
                                                            plot_fraction = 0.2, coverage_type = "line", 
                                                            rescale_introns = TRUE, heights = c(0.6,0.4)), .to = "plot")
plots_df_df = dplyr::mutate(plots_df, plot_title = paste(gene_name, snp_id, phenotype_id, sep = "_"))
plot_list = setNames(plots_df_df$plot, plots_df_df$plot_title)
savePlotList(plot_list, "results/coloc/plots/coverage/reviseAnnotations/")


#Visualise ensembl_87 annotations
ensembl_olaps = dplyr::select(gwas_olaps$ensembl_87, phenotype_id, snp_id, gene_name) %>% unique()
ensembl_hits = dplyr::semi_join(purrr::map_df(trqtl_min_pvalues$Ensembl_87, identity), ensembl_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::left_join(ensembl_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::select(group_id, phenotype_id, snp_id, gene_name) %>%
  unique()

#Make all coverage plots
plots_df = purrr::by_row(ensembl_hits, ~makeQTLCoveragePlot(.,str1_df, str2_df, vcf_file$genotypes,
                                                            exons,
                                                            gene_metadata = gene_meta_list$Ensembl_87,
                                                            plot_fraction = 0.2, coverage_type = "line", 
                                                            rescale_introns = TRUE, heights = c(0.6,0.4)), .to = "plot")
plots_df_df = dplyr::mutate(plots_df, plot_title = paste(gene_name, snp_id, phenotype_id, sep = "_"))
plot_list = setNames(plots_df_df$plot, plots_df_df$plot_title)
savePlotList(plot_list, "results/coloc/plots/coverage/Ensembl_87/")


