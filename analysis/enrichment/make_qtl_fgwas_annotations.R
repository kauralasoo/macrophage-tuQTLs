library("dplyr")
library("readr")
library("tximport")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import variant information
GRCh37_information = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")
var_info = dplyr::select(GRCh37_information, chr, pos, snp_id, MAF, type)

#Import QTL p-values
qtltools_columns = c("phenotype_id","pheno_chr","pheno_start", "pheno_end",
                    "strand","n_snps", "distance", "snp_id", "snp_chr",
                    "snp_start", "snp_end", "p_nominal","beta", "is_lead")
qtltools = readr::read_delim("processed/salmonella/fgwas/min_pvalues/naive.min_pvalues.txt.gz", delim = "\t", col_names = qtltools_columns) %>%
  dplyr::filter(p_nominal < 1e-4)
is_qtl_df = dplyr::transmute(qtltools, snp_id, featureCounts = 1)

#Add to var_info
new_annotations = dplyr::left_join(var_info, is_qtl_df, by = "snp_id") %>%
  dplyr::mutate(featureCounts = ifelse(is.na(featureCounts), 0, featureCounts))

#Save results
gz1 <- gzfile("processed/annotations/fgwas/fgwas_qtl_annotations.txt.gz", "w")
write.table(new_annotations, gz1, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz1)