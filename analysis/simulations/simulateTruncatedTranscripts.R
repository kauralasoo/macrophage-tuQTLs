library("dplyr")
library("BSgenome")
library("devtools")
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
  dplyr::left_join(truncated_transcripts, by = c("truncated_tx" = "ensembl_transcript_id")) %>%
  dplyr::mutate(truncation = case_when(
    cds_start_NF == 1 & cds_end_NF == 0 ~ "start",
    cds_start_NF == 0 & cds_end_NF == 1 ~ "end",
    cds_start_NF == 1 & cds_end_NF == 1 ~ "both"
  ))


#Calculate sequence differences in basepairs
findAllDiffs <- function(tx1, tx2, exons){
  print(paste(tx1, tx2))
  diff = txrevise::indentifyAddedRemovedRegions(tx1, tx2, exons) %>%
    calculateBasepairDifference()
}

#Find all differences between the two transcripts
tx1_list = as.list(truncated_pairs$full_tx)
tx2_list = as.list(truncated_pairs$truncated_tx)
all_differences = purrr::map2(tx1_list, tx2_list, ~findAllDiffs(.x, .y, exons)) %>% purrr::map_df(identity)

#Merge results
merged_diffs = dplyr::left_join(truncated_pairs, all_differences, by = c("full_tx" = "tx1_id")) %>% tbl_df()
unique_tx_ids = unique(c(merged_diffs$full_tx, merged_diffs$truncated_tx))

#Extract metadata for all transcripts
tx_meta = dplyr::filter(transcript_data, ensembl_transcript_id %in% unique_tx_ids) %>%
  txrevise::filterTranscriptMetadata()
saveRDS(tx_meta, "results/simulations/transcript_meta.rds")
tx_meta = readRDS("results/simulations/transcript_meta.rds")
tx_exons = exons[tx_meta$ensembl_transcript_id]
tx_cdss = cdss[intersect(tx_meta$ensembl_transcript_id, names(cdss))]

#Extend transcripts and construct events
extendTruncatedTx <- function(gene_id, tx_meta, exons, cdss){
  print(gene_id)
  
  #Extract gene data
  gene_data = txrevise::extractGeneData(gene_id, tx_meta, exons, cdss)
  
  #Extend transcripts
  gene_extended_tx = txrevise::extendTranscriptsPerGene(gene_data$metadata, gene_data$exons, gene_data$cdss)
  gene_data_ext = txrevise::replaceExtendedTranscripts(gene_data, gene_extended_tx)
  
  #Construct alt events
  alt_events = txrevise::constructAlternativeEvents(gene_data_ext$exons, gene_id)
  
  #Return results
  return(list(extended_tx = gene_data_ext, alt_events = alt_events))
}

#Apply to all genes
gene_ids = unique(tx_meta$ensembl_gene_id)
gene_ids_list = seqUtils::idVectorToList(gene_ids)
alt_events = purrr::map(gene_ids_list, ~extendTruncatedTx(., tx_meta, tx_exons, tx_cdss))
saveRDS(alt_events, "results/simulations/extended_tx_and_events.rds")
alt_events = readRDS("results/simulations/extended_tx_and_events.rds")

#Extract extended transcripts
new_exons = purrr::map(alt_events, ~as.list(.$extended_tx$exons)) %>% purrr::flatten()
new_exons = new_exons[names(tx_exons)]

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

old_exons_sorted = purrr::map(as.list(tx_exons), sortGrangesByStrand)
new_exons_sorted = purrr::map(new_exons, sortGrangesByStrand)

#Extract sequences
old_sequences = BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38, GRangesList(old_exons_sorted))
new_sequences = BSgenome::getSeq(BSgenome.Hsapiens.NCBI.GRCh38, GRangesList(new_exons_sorted))

#Concat exons
old_fastas = DNAStringSet(lapply(old_sequences, unlist))[tx_meta$ensembl_transcript_id]
new_fastas = DNAStringSet(lapply(new_sequences, unlist))[tx_meta$ensembl_transcript_id]

#Write transcripts to disk
writeXStringSet(old_fastas, 'results/simulations/original_transcripts.fa')
writeXStringSet(new_fastas, 'results/simulations/extended_transcripts.fa')


#Calculate effect sizes for tuQTLs
tx_meta = readRDS("results/simulations/transcript_meta.rds")
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
fasta = readDNAStringSet("results/simulations/original_transcripts.fa")
readspertx = round(50 * width(fasta) / 100)

simulate_experiment('results/simulations/original_transcripts.fa', reads_per_transcript=readspertx, 
                    num_reps=c(1,1), fold_changes=fold_changes[,1:2],
                    outdir='results/simulations/original_transcripts', gzip=TRUE, strand_specific = TRUE) 

#Simulate reads from the extended transcripts
fasta = readDNAStringSet("results/simulations/extended_transcripts.fa")
readspertx = round(50 * width(fasta) / 100)

simulate_experiment('results/simulations/extended_transcripts.fa', reads_per_transcript=readspertx, 
                    num_reps=rep(1,86), fold_changes=fold_changes,
                    outdir='results/simulations/extended_transcripts', gzip=TRUE, strand_specific = TRUE) 




#Construct alternative events
#Extend transcripts and construct events
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
alt_events = purrr::map(gene_ids_list, ~constructEvents(., tx_meta, tx_exons, tx_cdss))

#Flatten
alt_events = purrr::flatten(alt_events) %>% flattenAlternativeEvents()

#Construct event metadata
event_metadata = txrevise::constructEventMetadata(names(alt_events))

#Make annotations
annotations = txrevise::transcriptsToAnnotations(alt_events, event_metadata)
rtracklayer::export.gff3(annotations, "results/simulations/txrevise_annotatons.gff3")





