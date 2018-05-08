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


#Sort by chromosome and position
sorted_annotations = dplyr::arrange(new_annotations, chr, pos)


#Import GWAS summary stats
ra_gwas = readr::read_delim("RA_fgwas_input.txt.gz", delim = " ", col_types = "cciidid")

joint_annotations = dplyr::left_join(ra_gwas, sorted_annotations, by = c("CHR" = "chr", "POS" = "pos")) %>%
  dplyr::filter(!is.na(MAF))

filtered_annot = dplyr::select(joint_annotations, SNPID, CHR, POS, MAF, Z, N, SE, everything()) %>%
  dplyr::select(-type, -snp_id, -F) %>%
  dplyr::rename(F = MAF)

keep_unique_pos = dplyr::group_by(filtered_annot, CHR, POS) %>% dplyr::filter(row_number() == 1) %>% dplyr::ungroup()

#Save results
gz1 <- gzfile("RA_fgwas_annotated.txt.gz", "w")
write.table(keep_unique_pos, gz1, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz1)
