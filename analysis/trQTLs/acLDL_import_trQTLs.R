library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Import p-values
ensembl_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/Ensembl_87/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/Ensembl_87/AcLDL.permuted.txt.gz"))
revised_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/reviseAnnotations/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/reviseAnnotations/AcLDL.permuted.txt.gz"))
leafcutter_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/leafcutter/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/leafcutter/AcLDL.permuted.txt.gz"))
tpm_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/tpm/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/tpm/AcLDL.permuted.txt.gz"))
featureCounts_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/featureCounts/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/featureCounts/AcLDL.permuted.txt.gz"))
txrevise_promoters_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/txrevise_promoters/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/txrevise_promoters/AcLDL.permuted.txt.gz"))
txrevise_contained_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/txrevise_contained/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/txrevise_contained/AcLDL.permuted.txt.gz"))
txrevise_ends_pvalues = list(
  Ctrl = importQTLtoolsTable("processed/acLDL/qtltools/output/txrevise_ends/Ctrl.permuted.txt.gz"),
  AcLDL = importQTLtoolsTable("processed/acLDL/qtltools/output/txrevise_ends/AcLDL.permuted.txt.gz"))

#Put all results into a list
trQTL_min_pvalue_list = list(Ensembl_87 = ensembl_pvalues, reviseAnnotations = revised_pvalues, 
                             leafcutter = leafcutter_pvalues, tpm = tpm_pvalues,
                             featureCounts = featureCounts_pvalues, txrevise_promoters = txrevise_promoters_pvalues,
                             txrevise_contained = txrevise_contained_pvalues,
                             txrevise_ends = txrevise_ends_pvalues)
saveRDS(trQTL_min_pvalue_list, "results/trQTLs/acLDL_trQTL_min_pvalues.rds")


#Make qq-plots
qq_df = dplyr::mutate(ensembl_pvalues$Ctrl, p_eigen = p_beta) %>% dplyr::arrange(p_eigen) %>% addExpectedPvalue()
ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_eigen,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")