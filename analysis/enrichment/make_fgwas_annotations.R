library("dplyr")
library("readr")
library("rtracklayer")
library("tximport")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import variant information
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
variant_ranges = GRanges(seqnames = variant_information$chr, ranges = IRanges(start = variant_information$pos, end = variant_information$pos), strand = "*")
elementMetadata(variant_ranges) = as.data.frame(dplyr::select(variant_information, snp_id, MAF, type))

#Import ATAC peaks
atac_peaks = rtracklayer::import.gff3("processed/annotations/reference_QTLs/ATAC_consensus_peaks.gff3")
olaps = findOverlaps(variant_ranges, atac_peaks, ignore.strand=FALSE)
variant_ranges$atac_peak = 0
variant_ranges[queryHits(olaps)]$atac_peak = 1

#Export fgwas annotations
df = as.data.frame(variant_ranges) %>% 
  dplyr::select(-width, -strand, -end) %>%
  dplyr::rename(chr = seqnames, pos = start)

gz1 <- gzfile("processed/annotations/fgwas/fgwas_annotations.txt.gz", "w")
write.table(df, gz1)
close(gz1)
