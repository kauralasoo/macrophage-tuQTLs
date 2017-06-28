library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Import p-values
ensembl_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/IFNg_SL1344.permuted.txt.gz"),
  NI = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/NI.permuted.txt.gz"),
  NS = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/NS.permuted.txt.gz"),
  NIS = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/NIS.permuted.txt.gz"))
revised_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/IFNg_SL1344.permuted.txt.gz"),
  NI = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/NI.permuted.txt.gz"),
  NS = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/NS.permuted.txt.gz"),
  NIS = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/NIS.permuted.txt.gz"))
leafcutter_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/IFNg_SL1344.permuted.txt.gz"),
  NI = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/NI.permuted.txt.gz"),
  NS = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/NS.permuted.txt.gz"),
  NIS = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/NIS.permuted.txt.gz"))
tpm_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/IFNg_SL1344.permuted.txt.gz"),
  NI = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/NI.permuted.txt.gz"),
  NS = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/NS.permuted.txt.gz"),
  NIS = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/NIS.permuted.txt.gz"))
featureCounts_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/IFNg_SL1344.permuted.txt.gz"),
  NI = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/NI.permuted.txt.gz"),
  NS = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/NS.permuted.txt.gz"),
  NIS = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/NIS.permuted.txt.gz"))

#Put all results into a list
trQTL_min_pvalue_list = list(Ensembl_87 = ensembl_pvalues, reviseAnnotations = revised_pvalues, leafcutter = leafcutter_pvalues, tpm = tpm_pvalues,
                             featureCounts = featureCounts_pvalues)
saveRDS(trQTL_min_pvalue_list, "results/trQTLs/salmonella_trQTL_min_pvalues.rds")


#Make qq-plots
qq_df = dplyr::mutate(ensembl_pvalues$naive, p_eigen = p_beta) %>% dplyr::arrange(p_eigen) %>% addExpectedPvalue()
ggplot(qq_df, aes(x = -log(p_expected,10), y = -log(p_eigen,10))) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "black") + 
  theme_light() + 
  xlab("-log10 exptected p-value") + 
  ylab("-log10 observed p-value")

