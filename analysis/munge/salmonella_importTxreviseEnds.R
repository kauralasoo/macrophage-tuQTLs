library("dplyr")
library("readr")
library("tximport")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import Ensembl quant results for references
ensembl_quants = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")
sample_names = colnames(ensembl_quants)

#Iterate over annotations
annotations = c("txrevise.grp_1_ends","txrevise.grp_2_ends")
annotation_list = idVectorToList(annotations)

#Make lists of file names
file_names = purrr::map(annotation_list, ~setNames(file.path("processed/salmonella/salmon/",., sample_names, "quant.sf"), sample_names))

#Import all transcript abundances
tx_abundances = list('txrevise.grp_1_ends' = readRDS("processed/salmonella/matrices/txrevise.grp_1_ends.salmon_txrevise.rds"),
                     'txrevise.grp_2_ends' = readRDS("processed/salmonella/matrices/txrevise.grp_2_ends.salmon_txrevise.rds"))
tx_transposed = purrr::transpose(tx_abundances)

#Extract abundances, counts and lengths
abundances = purrr::reduce(tx_transposed$abundance, rbind)
counts = purrr::reduce(tx_transposed$counts, rbind)
lengths = purrr::reduce(tx_transposed$counts, rbind)

#Extract gene ids from event names
gene_meta = tbl_df2(rowData(ensembl_quants)) %>%
  dplyr::select(gene_id, gene_name, chr, strand, start, end) %>%
  unique()
gene_names = data_frame(transcript_id = rownames(abundances)) %>% 
  tidyr::separate(transcript_id, c("ensembl_gene_id", "group", "position", "ensembl_transcript_id"), "\\.", remove = FALSE) %>% 
  dplyr::mutate(gene_id = paste(ensembl_gene_id, position, sep = ".")) %>%
  dplyr::mutate(group_id = paste(ensembl_gene_id, group, position, sep = ".")) %>%
  dplyr::filter(ensembl_gene_id %in% gene_meta$gene_id) %>%
  dplyr::left_join(gene_meta, gene_meta, by = c("ensembl_gene_id" = "gene_id")) %>%
  dplyr::mutate(strand = ifelse(strand == 1, "+","-")) %>%
  dplyr::select(-group, -position) %>%
  as.data.frame()
rownames(gene_names) = gene_names$transcript_id

#Filter quants
abundances_filtered = abundances[gene_names$transcript_id,]
counts_filtered = counts[gene_names$transcript_id,]
lengths_filtered = lengths[gene_names$transcript_id,]

#Calculate abundance ratios
gene_name_map = dplyr::transmute(gene_names, gene_id = group_id, transcript_id)
abundance_ratios = calculateTranscriptRatios(abundances_filtered, gene_name_map)

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts_filtered, tpms = abundances_filtered, 
                relLengths = lengths_filtered, tpm_ratios = abundance_ratios), 
  colData = colData(ensembl_quants), 
  rowData = gene_names)

saveRDS(se, "results/SummarizedExperiments/salmonella_salmon_txrevise_ends.rds")



