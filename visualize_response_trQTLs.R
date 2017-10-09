library("dplyr")
library("devtools")
library("rtracklayer")
library("wiggleplotr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

#Import summary infromation for response trQTLs
response_trQTLs = readRDS("results/trQTLs/variance_explained/signifcant_interactions.rds") %>%
  dplyr::filter(quant != "featureCounts", condition != "AcLDL") 

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
revised_groups = dplyr::transmute(revised_gene_meta, phenotype_id = transcript_id, group_id)

##### Visualize IFNg+SL1344 specific trQTLs
#Construct dfs with paths to bigwig files
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

#Extract QTL hits
qtl_hits = dplyr::filter(response_trQTLs, condition == "IFNg_SL1344", quant == "reviseAnnotations", !all_technical) %>% 
  dplyr::left_join(revised_groups, by = "phenotype_id") %>%
  arrange(p_fdr)
hits_df = dplyr::select(qtl_hits, phenotype_id, snp_id, group_id, gene_name, condition, quant) %>%
  dplyr::mutate(plot_title = paste(condition, gene_name, snp_id, phenotype_id, sep = "_"))

#Make plots
plots_df = purrrlyr::by_row(hits_df, ~makeQTLCoveragePlot(., str1_df, str2_df, vcf_file$genotypes,
                    revised_granges, 
                    gene_metadata = revised_gene_meta,
                    plot_fraction = 0.2, coverage_type = "line", 
                    rescale_introns = TRUE, heights = c(0.6,0.4)))
plots_list = plots_df$.out
names(plots_list) = plots_df$plot_title
savePlotList(plots_list, "processed/salmonella/coverage_plots/reviseAnnotations/", suffix = ".pdf", width = 10, height = 10)


### TODO: make boxplots for those genes as well
#Visualize a couple of examples
SPOPL_data = constructQtlPlotDataFrame("ENSG00000138801.grp_1.contained.ENST00000512641", "rs4141322", 
                                       assays(se_revised)$tpm_ratios, 
                                       vcf_file$genotypes, 
                                       colData(se_revised) %>% tbl_df2(), 
                                       dplyr::mutate(revised_gene_meta, gene_id = transcript_id)) %>%
  dplyr::mutate(genotype_text = genotype_value)
plotQtlCol(SPOPL_data)