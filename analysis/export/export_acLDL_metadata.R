library("dplyr")
library("SummarizedExperiment")


#Import sample names
sample_names = read.table("analysis/data/sample_lists/acLDL_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Construct a design matrix from the sample names
design_matrix = constructDesignMatrix_acLDL(sample_names) %>%
  dplyr::filter(!(donor %in% c("mijn"))) %>%
  dplyr::arrange(donor, condition) %>%
  dplyr::filter(!(condition %in% c("LAL", "LAL_AcLDL"))) #Remove LAL samples

#Load sample metadata
sample_metadata = readRDS("analysis/data/covariates/compiled_acLDL_metadata.rds")
sample_meta = dplyr::left_join(design_matrix, sample_metadata, by = c("donor")) %>%
  dplyr::as.tbl() %>%
  dplyr::select(sample_id, line_id, replicate, condition_name, acLDL_date, ng_ul_mean, rna_extraction, rna_submit) %>%
  dplyr::mutate(chemistry = "V4", rna_auto = "FALSE")

#Import biosamples accessions
open = read.table("analysis/data/sample_lists/acLDL/biosamples/acLDL.open_access.txt", 
                  sep = ",", comment.char = "", stringsAsFactors = FALSE, header = TRUE) %>% tbl_df()
closed = read.table("analysis/data/sample_lists/acLDL/biosamples/acLDL.managed_access.txt", 
                  sep = ",", comment.char = "", stringsAsFactors = FALSE, header = TRUE) %>% tbl_df()
biosamples = dplyr::bind_rows(open,closed) %>% dplyr::transmute(sample_id, biosamples_accession = accession)
sample_meta = dplyr::left_join(sample_meta, biosamples, by = "sample_id")

write.table(sample_meta, "results/figures/tables/acLDL_sample_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)
