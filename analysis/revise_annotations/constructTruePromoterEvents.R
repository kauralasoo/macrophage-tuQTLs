library("purrr")
library("dplyr")
library("devtools")
library("rtracklayer")
library("GenomicFeatures")
library("optparse")
load_all("../txrevise/")
load_all("../seqUtils/")

#Read command-line options
option_list <- list(
  make_option(c("-t", "--transcripts"), type="character", default=NULL,
              help="Path to revised transcript annotations object.", metavar = "path"),
  make_option(c("-b", "--batch"), type="character", default = "NULL", 
              help = "Batch id.", metavar = "path"),
  make_option(c("-o", "--output"), type="character", default = "NULL", 
              help = "Path to the output folder.", metavar = "path"),
  make_option(c("-p", "--position"), type="character", default = "NULL", 
              help = "Make either promoter events ('start') or 3' end events ('end').", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Test options
#opt = list(t = "results/annotations/reviseAnnotations.GRangesList.rds",  b = "3 10000", o = "results/txrevise", p = "end")

#Import revised transcript annotations
revised_granges = readRDS(opt$t)

#Make gene metadata filter
if(opt$p == "start"){
  meta_filter = "upstresam"
} else if(opt$p == "end"){
  meta_filter = "downstream"
} else{
  error("Position (-p) has to be either 'start' or 'end'.")
}

#Construct gene metadata
events = txrevise::constructEventMetadata(names(revised_granges)) %>%
  dplyr::filter(event_type == meta_filter)
all_gene_ids = unique(events$gene_id)

#### Split genes into batches ####
batch_vector = as.integer(unlist(strsplit(opt$b, split = " ")))
gene_ids = unique(all_gene_ids)
batch_size = ceiling(length(gene_ids)/batch_vector[2])
batches = splitIntoBatches(length(gene_ids), batch_size)
selection = batches == batch_vector[1]
gene_ids = gene_ids[selection]

#Set up output file names
batch_id = paste(batch_vector, collapse = "_")
out_file = file.path(opt$o, paste0("txrevise_",opt$p,".batch_",batch_id, ".gff3"))
error_file = file.path(opt$o, paste0("failed_genes_",opt$p,".batch_",batch_id, ".txt"))

#Make a safe version of the function
safe_promoters = purrr::safely(fillMissingInternalExons)

#Construct new promotoers
alt_promoters = list()
failed_genes = c()
if(length(gene_ids) > 0){
  for(i in seq_along(gene_ids)){
    transcript_meta = dplyr::filter(events, gene_id == gene_ids[[i]])
    transcripts = revised_granges[transcript_meta$transcript_id]
    new_transcripts = safe_promoters(transcripts, type = opt$p)
    if(is.null(new_transcripts$error)){
      alt_promoters[[i]] = new_transcripts$result
    }else{
      alt_promoters[[i]] = transcripts
      failed_genes = c(failed_genes, gene_ids[[i]])
    }
  }
  
  #Flatten the results
  alt_promoters = purrr::reduce(alt_promoters, c)
  
  #Construct event metadata
  event_metadata = txrevise::constructEventMetadata(names(alt_promoters))
  
  #Construct trancript annotations
  annotations = txrevise::transcriptsToAnnotations(alt_promoters, event_metadata)
  
  #Save files
  rtracklayer::export.gff3(annotations, out_file)
  write.table(failed_genes, error_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
} else {
  write.table(c(), out_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(c(), error_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
}


