library("devtools")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("../wiggleplotr")

#Import SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")

#Extract TPM matrix
tpm_matrix = assays(se_revised)$tpm_ratios
TPM_matrix = assays(se_revised)$tpms
sample_meta = colData(se_revised) %>% tbl_df2()
gene_meta = rowData(se_revised) %>% tbl_df2() %>%
  dplyr::mutate(gene_id = transcript_id)

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
GRCh38_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
GRCh37_information = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")

#Import QTLs
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")


#Make a boxplot
IL1RN
boxplot = constructQtlPlotDataFrame("ENSG00000136689.grp_1.upstream.ENST00000409052", "rs28648961", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText("rs28648961", GRCh38_information), by = "genotype_value") %>%
  plotQtlCol()


boxplot = constructQtlPlotDataFrame("ENSG00000185507.grp_1.contained.ENST00000528413", "rs7101726", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText("rs7101726", GRCh38_information), by = "genotype_value") %>%
  plotQtlCol()


C12orf43
boxplot = constructQtlPlotDataFrame("ENSG00000185507.grp_1.contained.ENST00000528413", "rs7101726", 
                                    tpm_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText("rs7101726", GRCh38_information), by = "genotype_value") %>%
  plotQtlCol()

#CD40
boxplot = constructQtlPlotDataFrame("ENSG00000101017.grp_1.upstream.ENST00000372276", "rs4239702", 
                                    TPM_matrix, vcf_file$genotypes, sample_meta, gene_meta) %>% 
  dplyr::left_join(constructGenotypeText("rs4239702", GRCh38_information), by = "genotype_value") 
dplyr::filter(boxplot, condition_name == "IFNg")$norm_exp

%>%
  plotQtlCol()
