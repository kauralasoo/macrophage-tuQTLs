library("dplyr")
library("devtools")
library("ggplot2")
load_all("../seqUtils/")

#Import p-values
ensembl_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/Ensembl_87/IFNg_SL1344.permuted.txt.gz"))
revised_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations/IFNg_SL1344.permuted.txt.gz"))
leafcutter_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/leafcutter/IFNg_SL1344.permuted.txt.gz"))
tpm_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/tpm/IFNg_SL1344.permuted.txt.gz"))
featureCounts_pvalues = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/featureCounts/IFNg_SL1344.permuted.txt.gz"))
revised_groupwise = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations_groupwise/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations_groupwise/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations_groupwise/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/reviseAnnotations_groupwise/IFNg_SL1344.permuted.txt.gz"))
txrevise_promoters = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_promoters/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_promoters/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_promoters/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_promoters/IFNg_SL1344.permuted.txt.gz"))
txrevise_contained = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_contained/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_contained/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_contained/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_contained/IFNg_SL1344.permuted.txt.gz"))
txrevise_ends = list(
  naive = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_ends/naive.permuted.txt.gz"),
  IFNg = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_ends/IFNg.permuted.txt.gz"),
  SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_ends/SL1344.permuted.txt.gz"),
  IFNg_SL1344 = importQTLtoolsTable("processed/salmonella/qtltools/output/txrevise_ends/IFNg_SL1344.permuted.txt.gz"))

#Put all results into a list
trQTL_min_pvalue_list = list(Ensembl_87 = ensembl_pvalues, reviseAnnotations = revised_pvalues, leafcutter = leafcutter_pvalues, tpm = tpm_pvalues,
                             featureCounts = featureCounts_pvalues, reviseAnnotations_groupwise = revised_groupwise, txrevise_promoters = txrevise_promoters, 
                             txrevise_contained = txrevise_contained,
                             txrevise_ends = txrevise_ends)
saveRDS(trQTL_min_pvalue_list, "results/trQTLs/salmonella_trQTL_min_pvalues.rds")






