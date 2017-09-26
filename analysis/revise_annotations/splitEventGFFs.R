library("rtracklayer")
library("dplyr")
library("devtools")
load_all("../reviseAnnotations/")
load_all("../seqUtils/")
library("optparse")

#Read command-line options
option_list <- list(
  make_option(c("-a", "--gff1"), type="character", default=NULL,
              help="Path to the first gff file.", metavar = "path"),
  make_option(c("-b", "--gff2"), type="character", default=NULL,
              help="Path to the second gff file.", metavar = "path"),
  make_option(c("-o", "--output"), type="character", default = "NULL", 
              help = "Path to the output folder.", metavar = "path")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Imort the full GFF file
gff_grp1 = rtracklayer::import.gff3(opt$gff1)
gff_grp2 = rtracklayer::import.gff3(opt$gff2)

#Construct metadata
tx_records = dplyr::bind_rows(as.data.frame(gff_grp1) %>% tbl_df(), as.data.frame(gff_grp2) %>% tbl_df()) %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::transmute(transcript_id = ID, transcript_start = start, transcript_end = end, 
                   transcript_length = width, strand = as.character(strand), chr = as.character(seqnames))
event_metadata = reviseAnnotations::constructEventMetadata(tx_records$transcript_id)
tx_metadata = dplyr::left_join(tx_records, event_metadata, by = "transcript_id")
write.table(tx_metadata, file.path(opt$output, "reviseAnnotations.transcript_metadata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#split events by position
grp_1_upstream_gff = gff_grp1[grep("upstream", gff_grp1$gene_id),]
grp_1_downstream_gff = gff_grp1[grep("downstream", gff_grp1$gene_id),]
grp_1_contained_gff = gff_grp1[grep("contained", gff_grp1$gene_id),]

#Export different events
rtracklayer::export.gff3(grp_1_upstream_gff, file.path(opt$output, "reviseAnnotations.grp_1_upstream.gff3"))
rtracklayer::export.gff3(grp_1_downstream_gff, file.path(opt$output, "reviseAnnotations.grp_1_downstream.gff3"))
rtracklayer::export.gff3(grp_1_contained_gff, file.path(opt$output, "reviseAnnotations.grp_1_contained.gff3"))

#split events by position
grp_2_upstream_gff = gff_grp2[grep("upstream", gff_grp2$gene_id),]
grp_2_downstream_gff = gff_grp2[grep("downstream", gff_grp2$gene_id),]
grp_2_contained_gff = gff_grp2[grep("contained", gff_grp2$gene_id),]

#Export different events
rtracklayer::export.gff3(grp_2_upstream_gff, file.path(opt$output, "reviseAnnotations.grp_2_upstream.gff3"))
rtracklayer::export.gff3(grp_2_downstream_gff, file.path(opt$output, "reviseAnnotations.grp_2_downstream.gff3"))
rtracklayer::export.gff3(grp_2_contained_gff, file.path(opt$output, "reviseAnnotations.grp_2_contained.gff3"))


