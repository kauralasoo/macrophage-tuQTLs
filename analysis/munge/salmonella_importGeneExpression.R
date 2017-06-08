library("dplyr")
library("readr")
library("tximport")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import sample names
sample_names = read.table("analysis/data/sample_lists/SL1344_names_all.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE)[,1]

#Construct a design matrix from the sample names
design_matrix = constructDesignMatrix_SL1344(sample_names) %>% #Construct a design matrix from the sample names
  dplyr::filter(!(donor == "fpdj")) %>% tbl_df() %>% #Remove both fpdj samples (same as nibo)
  dplyr::filter(!(donor == "fpdl" & replicate == 2)) %>% #Remove second fpdl sample (ffdp)
  dplyr::filter(!(donor == "ougl" & replicate == 2)) %>% #Remove second ougl sample (dium)
  dplyr::filter(!(donor == "mijn")) %>% #Remove mijn (wrong line from CGAP)
  dplyr::filter(!(donor == "qaqx")) %>% #Remove qaqx (wrong line from CGAP) 
  dplyr::mutate(replicate = ifelse(donor == "babk",2,replicate)) %>% #Change babk replicate to two
  dplyr::arrange(donor, condition) %>%
  as.data.frame()

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds")) %>%
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, gene_name = external_gene_name, chr = chromosome_name)
filtered_transcscript_data = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.filtered.rds")

#Create a vector of file names
file_names = file.path("processed/salmonella/salmon/Ensembl_87/", design_matrix$sample_id, "quant.sf")
names(file_names) = design_matrix$sample_id

#Convert transcript data to suitable format for tximport
tx2gene = dplyr::select(transcript_data, gene_id, transcript_id, transcript_version) %>% 
  dplyr::mutate(TXNAME = paste(transcript_id, transcript_version, sep = ".")) %>% 
  dplyr::transmute(TXNAME,gene_id, transcript_id)

#Import gene-level abundances
gene_abundances = tximport(file_names, type = "salmon", tx2gene = tx2gene[,1:2], importer = read_tsv, dropInfReps = TRUE)

#Prepare sample metadata
sample_metadata = readRDS("analysis/data/covariates/compiled_salmonella_metadata.rds")
sample_meta = dplyr::left_join(design_matrix, sample_metadata, by = c("donor", "replicate")) %>% as.data.frame()
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
saveRDS(se, "results/SummarizedExperiments/salmonella_salmon_gene_abundances.rds")

