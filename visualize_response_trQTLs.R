library("dplyr")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

#Import summary infromation for response trQTLs
salmonella_varExp = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")

a = dplyr::filter(salmonella_varExp, condition == "IFNg_SL1344", quant == "reviseAnnotations") %>%
  dplyr::arrange(desc(interaction_fraction))

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")

#Import gene exrression data
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")

#Import Ensembl transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

#Import revisedAnnotations Granges
revised_granges = readRDS("results/annotations/reviseAnnotations.GRangesList.rds")
leafcutter_granges = readRDS("results/annotations/salmonella_leafcutter.GRangesList.rds")

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
revised_gene_meta = rowData(se_revised) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = group_id)

#Set up sample coverage df
str1_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                       colData(se_revised) %>% tbl_df2(), 
                                       "/Volumes/Ajamasin/bigwig/salmonella/", 
                                       bigWig_suffix = ".str1.bw",
                                       condition_name_levels = c("naive","IFNg", "SL1344", "IFNg_SL1344")) %>%
  dplyr::filter(condition_name %in% c("naive", "IFNg_SL1344"))
str2_df = wiggleplotrConstructMetadata(assays(se_featureCounts)$counts, 
                                       colData(se_revised) %>% tbl_df2(), 
                                       "/Volumes/Ajamasin/bigwig/salmonella/", 
                                       bigWig_suffix = ".str2.bw",
                                       condition_name_levels = c("naive","IFNg", "SL1344", "IFNg_SL1344")) %>%
  dplyr::filter(condition_name %in% c("naive", "IFNg_SL1344"))

#QTL hits
a = dplyr::filter(salmonella_varExp, condition == "IFNg_SL1344", quant == "reviseAnnotations") %>%
  dplyr::arrange(desc(interaction_fraction)) %>%
  dplyr::left_join(dplyr::transmute(revised_gene_meta, phenotype_id = transcript_id, group_id))
hit = dplyr::select(a, phenotype_id, snp_id, group_id) %>% dplyr::filter(row_number() == 25)

makeQTLCoveragePlot(hit, str1_df, str2_df, vcf_file$genotypes,
                    revised_granges, 
                    gene_metadata = revised_gene_meta,
                    plot_fraction = 0.2, coverage_type = "line", 
                    rescale_introns = TRUE, heights = c(0.6,0.4))

#Import QTL mapping results
trqtl_min_pvalues = readRDS("results/acLDL/trQTLs/trQTL_min_pvalues.other_tx.rds")

#Import SummarizedExperiments for all phenotypes
se_ensembl = readRDS("results/acLDL/acLDL_salmon_ensembl.rds")
se_revised = readRDS("results/acLDL/acLDL_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/acLDL/acLDL_leafcutter_counts.rds")
se_list = list(ensembl_87 = se_ensembl, revisedAnnotation = se_revised, leafcutter = se_leafcutter)

#Extract sample metadata
sample_meta_list = purrr::map(se_list, ~colData(.) %>% tbl_df2())
gene_meta_list = purrr::map(se_list, ~rowData(.) %>% tbl_df2())


#Set up sample coverage df
str1_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_meta_list$revisedAnnotation, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str1.bw",
                                       condition_name_levels = c("Ctrl","AcLDL"))
str2_df = wiggleplotrConstructMetadata(acldl_list$counts, 
                                       sample_meta_list$revisedAnnotation, 
                                       "/Volumes/JetDrive/bigwig/", 
                                       bigWig_suffix = ".str2.bw",
                                       condition_name_levels = c("Ctrl","AcLDL"))



#### Visualise GWAS overlaps ####
gwas_olaps = readRDS("acLDL_figures/tables/GWAS_coloc_hits.rds")

#Visualise leafCutter junctions
leafcutter_olaps = dplyr::select(gwas_olaps$leafcutter, phenotype_id, snp_id, gene_name) %>% unique()
leafcutter_hits = dplyr::semi_join(purrr::map_df(trqtl_min_pvalues$leafcutter, identity), leafcutter_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::left_join(leafcutter_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::select(group_id,phenotype_id, other_phenotype_id, snp_id, gene_name) %>%
  unique()

#Make all coverage plots
plots_df = purrr::by_row(leafcutter_hits, ~makeQTLCoveragePlot(.,str1_df, str2_df, vcf_file$genotypes,
                                                               leafcutter_granges, 
                                                               gene_metadata = gene_meta_list$leafcutter,
                                                               plot_fraction = 0.2, coverage_type = "line", 
                                                               rescale_introns = TRUE, heights = c(0.6,0.4)), .to = "plot")
plots_df_df = dplyr::mutate(plots_df, plot_title = paste(gene_name, snp_id, phenotype_id, sep = "_"))
plot_list = setNames(plots_df_df$plot, plots_df_df$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/coverage/leafcutter/")


#Visualise revisedAnnotation
revised_olaps = dplyr::select(gwas_olaps$revisedAnnotation, phenotype_id, snp_id, gene_name) %>% unique()
revised_hits = dplyr::semi_join(purrr::map_df(trqtl_min_pvalues$revisedAnnotations, identity), revised_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::left_join(revised_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::select(group_id, phenotype_id, other_phenotype_id, snp_id, gene_name) %>%
  unique()

#Make all coverage plots
plots_df = purrr::by_row(revised_hits, ~makeQTLCoveragePlot(.,str1_df, str2_df, vcf_file$genotypes,
                                                            revised_granges, 
                                                            gene_metadata = gene_meta_list$revisedAnnotation,
                                                            plot_fraction = 0.2, coverage_type = "line", 
                                                            rescale_introns = TRUE, heights = c(0.6,0.4)), .to = "plot")
plots_df_df = dplyr::mutate(plots_df, plot_title = paste(gene_name, snp_id, phenotype_id, sep = "_"))
plot_list = setNames(plots_df_df$plot, plots_df_df$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/coverage/revisedAnnotation/")


#Visualise ensembl_87 annotations
ensembl_olaps = dplyr::select(gwas_olaps$ensembl_87, phenotype_id, snp_id, gene_name) %>% unique()
ensembl_hits = dplyr::semi_join(purrr::map_df(trqtl_min_pvalues$ensembl_87, identity), ensembl_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::left_join(ensembl_olaps, by = c("phenotype_id", "snp_id")) %>%
  dplyr::select(group_id, phenotype_id, other_phenotype_id, snp_id, gene_name) %>%
  unique()

#Make all coverage plots
plots_df = purrr::by_row(ensembl_hits, ~makeQTLCoveragePlot(.,str1_df, str2_df, vcf_file$genotypes,
                                                            exons,
                                                            gene_metadata = gene_meta_list$ensembl_87,
                                                            plot_fraction = 0.2, coverage_type = "line", 
                                                            rescale_introns = TRUE, heights = c(0.6,0.4)), .to = "plot")
plots_df_df = dplyr::mutate(plots_df, plot_title = paste(gene_name, snp_id, phenotype_id, sep = "_"))
plot_list = setNames(plots_df_df$plot, plots_df_df$plot_title)
savePlotList(plot_list, "processed/acLDL/coloc_plots/coverage/ensembl_87//")



#Visualise interaction resluts
interaction_list = readRDS("results/acLDL/trQTLs/trQTL_interaction_results.rds")
a = purrr::map_df(trqtl_min_pvalues$leafcutter, identity) %>% dplyr::filter(phenotype_id == "17:81192593:81192736:clu_27768") %>% head(n=1)

makeQTLCoveragePlot(a,str1_df, str2_df, vcf_file$genotypes,
                    leafcutter_granges, 
                    gene_metadata = gene_meta_list$leafcutter,
                    plot_fraction = 0.2, coverage_type = "line", 
                    rescale_introns = TRUE, heights = c(0.6,0.4))