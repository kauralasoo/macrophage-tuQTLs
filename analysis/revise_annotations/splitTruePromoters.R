library("rtracklayer")
library("dplyr")
library("devtools")
library("data.table")
load_all("../txrevise/")
load_all("../seqUtils/")

new_promoters = rtracklayer::import.gff3("processed/annotations/txrevise_promoters/merged/txrevise_promoters.gff3.gz")

tx_records = as.data.frame(new_promoters) %>% tbl_df() %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::transmute(transcript_id = ID, transcript_start = start, transcript_end = end, 
                   transcript_length = width, strand = as.character(strand), chr = as.character(seqnames))

#Expract groups
promoters_grp1 = new_promoters[new_promoters$gene_id %like% "grp_1"]
promoters_grp2 = new_promoters[new_promoters$gene_id %like% "grp_2"]

#Export
out_dir = "processed/annotations/gff/"
rtracklayer::export.gff3(promoters_grp1, file.path(out_dir, "txrevise.grp_1_promoters.gff3"))
rtracklayer::export.gff3(promoters_grp2, file.path(out_dir, "txrevise.grp_2_promoters.gff3"))


#Split alternative ends
new_promoters = rtracklayer::import.gff3("processed/annotations/txrevise_promoters/merged/txrevise_end.gff3.gz")

tx_records = as.data.frame(new_promoters) %>% tbl_df() %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::transmute(transcript_id = ID, transcript_start = start, transcript_end = end, 
                   transcript_length = width, strand = as.character(strand), chr = as.character(seqnames))

#Expract groups
promoters_grp1 = new_promoters[new_promoters$gene_id %like% "grp_1"]
promoters_grp2 = new_promoters[new_promoters$gene_id %like% "grp_2"]

#Export
out_dir = "processed/annotations/gff/"
rtracklayer::export.gff3(promoters_grp1, file.path(out_dir, "txrevise.grp_1_ends.gff3"))
rtracklayer::export.gff3(promoters_grp2, file.path(out_dir, "txrevise.grp_2_ends.gff3"))

