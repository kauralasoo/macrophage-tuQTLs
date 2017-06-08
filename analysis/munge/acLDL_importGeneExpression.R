library("dplyr")
library("readr")
library("tximport")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import sample names
sample_names = read.table("analysis/data/sample_lists/acLDL_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Construct a design matrix from the sample names
design_matrix = constructDesignMatrix_acLDL(sample_names) %>%
  dplyr::filter(!(donor %in% c("mijn", "xegx"))) %>% #Remove two samples (xeqx failed, mijn given to us by accident)
  dplyr::arrange(donor, condition) %>%
  dplyr::filter(!(condition %in% c("LAL", "LAL_AcLDL"))) #Remove LAL samples

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds")) %>%
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, gene_name = external_gene_name, chr = chromosome_name)
filtered_transcscript_data = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.filtered.rds")

#Create a vector of file names
file_names = file.path("processed/acLDL/salmon/Ensembl_87/", design_matrix$sample_id, "quant.sf")
names(file_names) = design_matrix$sample_id

#Convert transcript data to suitable format for tximport
tx2gene = dplyr::select(transcript_data, gene_id, transcript_id, transcript_version) %>% 
  dplyr::mutate(TXNAME = paste(transcript_id, transcript_version, sep = ".")) %>% 
  dplyr::transmute(TXNAME,gene_id, transcript_id)

#Import gene-level abundances
gene_abundances = tximport(file_names, type = "salmon", tx2gene = tx2gene[,1:2], importer = read_tsv, dropInfReps = TRUE)

#Prepare sample metadata
sample_metadata = readRDS("analysis/data/covariates/compiled_acLDL_metadata.rds")
sample_meta = dplyr::left_join(design_matrix, sample_metadata, by = c("donor")) %>% as.data.frame()
rownames(sample_meta) = sample_meta$sample_id

#Construct transcript metadata
gene_metadata = dplyr::filter(transcript_data, gene_id %in% rownames(gene_abundances$abundance)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>%
  dplyr::select(gene_id, gene_biotype, chr, start, end, gene_name, strand) %>% 
  unique() %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(gene_metadata) = gene_metadata$gene_id

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = gene_abundances$counts, tpms = gene_abundances$abundance, 
                relLengths = gene_abundances$length), 
  colData = sample_meta, 
  rowData = gene_metadata)
saveRDS(se, "results/SummarizedExperiments/acLDL_salmon_gene_abundances.rds")

