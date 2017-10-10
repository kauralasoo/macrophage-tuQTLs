library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("../seqUtils/")

#Parse command-line options
option_list <- list(
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("-g", "--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("-q", "--qtl"), type="character", default=NULL,
              help="Path to the QTL directory.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("-s", "--samplesizes"), type="character", default=NULL,
              help="Path to the tab-separated text file with condition names and sample sizes.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
#opt = list(g = "IBD", w = "2e5", p = "featureCounts", d = "/Volumes/JetDrive/datasets/Inflammatory_GWAS/", o = "results/acLDL/coloc/coloc_lists/",
#           q = "processed/salmonella/qtltools/output/", s = "analysis/data/sample_lists/salmonella_coloc_sample_sizes.txt")

#Extract parameters for CMD options
gwas_id = opt$g
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
qtl_dir = opt$q
outdir = opt$o
sample_size_path = opt$s

#Import variant information
GRCh38_variants = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
GRCh37_variants = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))

#Import sample sizes
sample_sizes = readr::read_tsv(sample_size_path, col_names = c("condition_name", "sample_size"), col_types = "cc")
sample_sizes_list = as.list(sample_sizes$sample_size)
names(sample_sizes_list) = sample_sizes$condition_name

#Construct a new QTL list 
phenotype_values = constructQtlListForColoc(phenotype, qtl_dir, sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, 
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#Test for coloc
coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes, 
                             ~colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x, 
                                                      gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
                                                      GRCh37_variants = GRCh37_variants,
                                                      GRCh38_variants = GRCh38_variants, 
                                                      N_qtl = .y, cis_dist = cis_window))

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead)
coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w, "txt", sep = "."))
write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging example
#colocMolecularQTLs(qtl_pairs[1,], qtl_summary_list$Ctrl, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants, GRCh38_variants, N_qtl = 2e5)
