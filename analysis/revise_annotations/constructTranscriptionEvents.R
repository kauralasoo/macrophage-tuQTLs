library("purrr")
library("dplyr")
library("devtools")
library("rtracklayer")
library("GenomicFeatures")
library("optparse")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")

#Read command-line options
option_list <- list(
  make_option(c("-t", "--txdb"), type="character", default=NULL,
              help="Path to the TxDb object.", metavar = "path"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Path to gene metadata table.", metavar = "path"),
  make_option(c("-b", "--batch"), type="character", default = "NULL", 
              help = "Batch id.", metavar = "path"),
  make_option(c("-o", "--output"), type="character", default = "NULL", 
              help = "Path to the output folder.", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Test options
#opt = list(t = "../../annotations/GRCh37/genes/Ensembl_90/TranscriptDb_GRCh37_90.db", 
#           m = "../../annotations/GRCh37/genes/Ensembl_90/Homo_sapiens.GRCh37.90.compiled_tx_metadata.txt.gz",
#           b = "2 10000",
#           o = "results/reviseAnnotation")

#Import transcript annotations
gene_metadata = readr::read_tsv(opt$m)
txdb = AnnotationDbi::loadDb(opt$t)
exons = GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

#Filter gene annotations
filtered_metadata = reviseAnnotations::filterTranscriptMetadata(gene_metadata)

#### Split genes into batches ####
batch_vector = as.integer(unlist(strsplit(opt$b, split = " ")))
gene_ids = unique(filtered_metadata$ensembl_gene_id)
batch_size = ceiling(length(gene_ids)/batch_vector[2])
batches = splitIntoBatches(length(gene_ids), batch_size)
selection = batches == batch_vector[1]
gene_ids = gene_ids[selection]

#Set up output file names
batch_id = paste(batch_vector, collapse = "_")
grp1_file = file.path(opt$o, paste0("reviseAnnotations.grp_1.batch_",batch_id, ".gff3"))
grp2_file = file.path(opt$o, paste0("reviseAnnotations.grp_2.batch_",batch_id, ".gff3"))
error_file = file.path(opt$o, paste0("failed_genes.batch_",batch_id, ".txt"))

#Only proceed with event construction if there are any genes in the list
if (length(gene_ids) > 0){
  #Construct events
  gene_ids_list = seqUtils::idVectorToList(gene_ids)
  
  #Construct alternative events and remove failed genes
  safe_construct = purrr::safely(constructAlternativeEventsWrapper)
  alt_events = purrr::map(gene_ids_list, ~safe_construct(., filtered_metadata, exons, cdss)$result)
  failed_genes = purrr::map_lgl(alt_events, is.null)
  alt_events = alt_events[!failed_genes] #Remove failed genes
  
  #Flatten
  alt_events = purrr::flatten(alt_events) %>% flattenAlternativeEvents()
  
  #Construct event metadata
  event_metadata = reviseAnnotations::constructEventMetadata(names(alt_events))
  
  #Separate the two groups
  grp1_events = dplyr::filter(event_metadata, grp_id == "grp_1")
  grp2_events = dplyr::filter(event_metadata, grp_id == "grp_2")
  
  #Construct transcript annotations
  grp1_annotations = reviseAnnotations::transcriptsToAnnotations(alt_events[grp1_events$transcript_id], grp1_events)
  grp2_annotations = reviseAnnotations::transcriptsToAnnotations(alt_events[grp2_events$transcript_id], grp2_events)
  
  #Make a list of failed genes
  failed_names = names(which(failed_genes))
  
  #Save output files to disk
  rtracklayer::export.gff3(grp1_annotations, grp1_file)
  rtracklayer::export.gff3(grp2_annotations, grp2_file)
  write.table(failed_names, error_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
} else {
  #Make empty output files
  write.table(c(), grp1_file, quote = FALSE)
  write.table(c(), grp2_file, quote = FALSE)
  write.table(c(), error_file, quote = FALSE)
}
