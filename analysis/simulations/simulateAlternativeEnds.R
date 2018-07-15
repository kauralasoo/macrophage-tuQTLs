library("dplyr")
library("polyester")
library("BSgenome")
library("devtools")
library("data.table")
library("GenomicRanges")
library("GenomicFeatures")
load_all("../txrevise/")
load_all("../seqUtils/")
library("BSgenome.Hsapiens.NCBI.GRCh38")

#Import transcript annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
exons = exonsBy(txdb, by = "tx", use.names=TRUE)
cdss = cdsBy(txdb, by = "tx", use.names=TRUE)

#Import QTLs
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")

#Import QTL pairs
QTL_pairs = readRDS("results/simulations/trQTL_pair_diffs.rds")

#Import transcript metadata
transcript_data = tbl_df(readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds"))
transcript_meta = dplyr::select(transcript_data, ensembl_transcript_id, cds_start_NF, cds_end_NF)

truncated_transcripts = dplyr::filter(transcript_meta, cds_start_NF == 1 | cds_end_NF == 1)

#Identify trQTL pairs with truncated transcripts
first_truncated = dplyr::semi_join(QTL_pairs, truncated_transcripts, by = c("tx1_id" = "ensembl_transcript_id"))
second_truncated = dplyr::semi_join(QTL_pairs, truncated_transcripts, by = c("tx2_id" = "ensembl_transcript_id"))
second_nonoverlap = dplyr::anti_join(second_truncated, first_truncated, by = "tx1_id")

#Make trancript pairs
truncated_pairs = dplyr::bind_rows(dplyr::transmute(first_truncated, full_tx = tx2_id, truncated_tx = tx1_id),
                                   dplyr::transmute(second_nonoverlap, full_tx = tx1_id, truncated_tx = tx2_id)) %>%
  as_tibble()

#Extract full-length transcripts
full_length_txs = exons[truncated_pairs$full_tx]
lengths = lapply(full_length_txs, length)
full_length_txs = full_length_txs[as.numeric(lengths) >= 5]
tx_pairs = data_frame(full_tx = names(full_length_txs)) %>% 
  dplyr::left_join(truncated_pairs)
unique_tx_ids = unique(c(tx_pairs$full_tx, tx_pairs$truncated_tx))

#Extract metadata for all transcripts
tx_meta = dplyr::filter(transcript_data, ensembl_transcript_id %in% unique_tx_ids) %>%
  txrevise::filterTranscriptMetadata()

#Remove 20% of the exons from both ends
removeExons <- function(tx, remove_fraction = 0.2, both = FALSE){
  tx = sort(tx)
  exon_count = length(tx)
  remove_count = floor(0.2*exon_count)
  
  if(both){
    result = tx[(remove_count+1):(length(tx)-remove_count)]
  } else {
    result = tx[(remove_count+1):length(tx)]
  }
  return(result)
}

#Remove both ends
truncated_both = purrr::map(as.list(full_length_txs), ~removeExons(., both = TRUE))
names(truncated_both) = tx_pairs$truncated_tx

#Remove first end
truncated_one = purrr::map(as.list(full_length_txs), ~removeExons(., both = FALSE))
names(truncated_one) = tx_pairs$truncated_tx

#Sort exons by strand
sortGrangesByStrand <- function(granges){
  tx_strand = as.character(strand(granges))[1]
  if(tx_strand == "-"){
    granges = sort(granges, decreasing = T)
  } else{
    granges = sort(granges, decreasing = F)
  }
  return(granges)
}

one_tx_sorted = c(as.list(full_length_txs), truncated_one) %>%
  purrr::map(sortGrangesByStrand)
both_tx_sorted = c(as.list(full_length_txs), truncated_both) %>%
  purrr::map(sortGrangesByStrand)

#Extract sequences
one_sequences = BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38, GRangesList(one_tx_sorted))
both_sequences = BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38, GRangesList(both_tx_sorted))

#Concat exons
one_fastas = DNAStringSet(lapply(one_sequences, unlist))[tx_meta$ensembl_transcript_id]
both_fastas = DNAStringSet(lapply(both_sequences, unlist))[tx_meta$ensembl_transcript_id]

#Write transcripts to disk
writeXStringSet(one_fastas, 'results/simulations/one_transcripts.fa')
writeXStringSet(both_fastas, 'results/simulations/both_transcripts.fa')


#Calculate effect sizes for tuQTLs
lead_variants = dplyr::filter(salmonella_qtls$Ensembl_87$naive, group_id %in% tx_meta$ensembl_gene_id) %>%
  dplyr::select(group_id, snp_id)
genotype_matrix = vcf_file$genotypes[lead_variants$snp_id,]
genotype_df = as_tibble(genotype_matrix) %>%
  dplyr::mutate(ensembl_gene_id = lead_variants$group_id) %>%
  dplyr::select(ensembl_gene_id, everything())

#Add effect size multiplier to tx_meta
set.seed(1)
effect_direction = dplyr::select(tx_meta, ensembl_gene_id, ensembl_transcript_id) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(effect_multiplier = c(1,-1)) %>% 
  dplyr::mutate(is_de = round(runif(1,0,1))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(effect = effect_multiplier*is_de) %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, effect)

#Make effect size matrix
fc_matrix = dplyr::left_join(effect_direction, genotype_df, by = "ensembl_gene_id") %>%
  dplyr::select(-ensembl_gene_id, -ensembl_transcript_id, -effect) %>%
  as.matrix()
row.names(fc_matrix) = effect_direction$ensembl_transcript_id
fc_matrix = effect_direction$effect*fc_matrix
fold_changes = 2^fc_matrix
fold_changes[is.na(fold_changes)] = 1

#Simulate reads from the original transcripts
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
fasta = readDNAStringSet("results/simulations/one_transcripts.fa")
readspertx = round(50 * width(fasta) / 100)

simulate_experiment('results/simulations/one_transcripts.fa', reads_per_transcript=readspertx, 
                    num_reps=rep(1,86), fold_changes=fold_changes,
                    outdir='results/simulations/one_transcripts', gzip=TRUE, strand_specific = TRUE) 

#Simulate reads from the extended transcripts
fasta = readDNAStringSet("results/simulations/both_transcripts.fa")
readspertx = round(50 * width(fasta) / 100)

simulate_experiment('results/simulations/both_transcripts.fa', reads_per_transcript=readspertx, 
                    num_reps=rep(1,86), fold_changes=fold_changes,
                    outdir='results/simulations/both_transcripts', gzip=TRUE, strand_specific = TRUE) 

#Construct alternative events
constructEvents <- function(gene_id, tx_meta, exons, cdss){
  print(gene_id)
  
  #Extract gene data
  gene_data = txrevise::extractGeneData(gene_id, tx_meta, exons, cdss)
  
  #Construct alt events
  alt_events = txrevise::constructAlternativeEvents(gene_data$exons, gene_id)
  
  #Return results
  return(alt_events)
}

#Apply to all genes
gene_ids = unique(tx_meta$ensembl_gene_id)
gene_ids_list = seqUtils::idVectorToList(gene_ids)
tx_exons = GRangesList(both_tx_sorted)
alt_events = purrr::map(gene_ids_list, ~constructEvents(., tx_meta, tx_exons, tx_exons))
saveRDS(alt_events, "results/simulations/one_both_alt_events.rds")

#Flatten
alt_events_flat = purrr::flatten(alt_events) %>% flattenAlternativeEvents()

#Construct event metadata
event_metadata = txrevise::constructEventMetadata(names(alt_events_flat))

#Make annotations
annotations = txrevise::transcriptsToAnnotations(alt_events_flat, event_metadata)
rtracklayer::export.gff3(annotations[annotations$gene_id %like% "upstream"], "results/simulations/txrevise_upstream_both.gff3")
rtracklayer::export.gff3(annotations[annotations$gene_id %like% "downstream"], "results/simulations/txrevise_downstream_both.gff3")



