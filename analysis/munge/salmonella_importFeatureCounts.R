library("dplyr")
library("readr")
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

#Import featureCounts
count_matrix = loadCounts("processed/salmonella/featureCounts/", design_matrix$sample_id, sub_dir = FALSE, counts_suffix = ".featureCounts.txt")

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.transcript_data.rds")) %>%
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, gene_name = external_gene_name, chr = chromosome_name)

#Filter transcript metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
filtered_tx_data = dplyr::filter(transcript_data, gene_biotype %in% valid_gene_biotypes, chr %in% valid_chromosomes)

#Extract length
length_df = dplyr::select(count_matrix, gene_id, length)

#Construct transcript metadata
gene_metadata = filtered_tx_data %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>%
  dplyr::select(gene_id, gene_biotype, chr, start, end, gene_name, strand, percentage_gc_content) %>% 
  dplyr::left_join(length_df, by = "gene_id") %>%
  unique() %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(gene_metadata) = gene_metadata$gene_id

#Process counts
filtered_data = dplyr::filter(count_matrix, gene_id %in% gene_metadata$gene_id)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id
counts = counts[gene_metadata$gene_id,] #Reoder counts

#CQN normalize counts
cqn_matrix = calculateCQN(counts, gene_metadata)

#Prepare sample metadata
sample_metadata = readRDS("analysis/data/covariates/compiled_salmonella_metadata.rds")
sample_meta = dplyr::left_join(design_matrix, sample_metadata, by = c("donor", "replicate")) %>% as.data.frame()
rownames(sample_meta) = sample_meta$sample_id

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts, cqn = cqn_matrix), 
  colData = sample_meta, 
  rowData = gene_metadata)
saveRDS(se, "results/SummarizedExperiments/salmonella_featureCounts.rds")




