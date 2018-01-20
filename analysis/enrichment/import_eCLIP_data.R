library("dplyr")
library("tidyr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")

valid_chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chrY", "chr19", "chr22", "chr21")

#Import metadata
metadata = read.table("processed/annotations/fgwas/raw/eCLIP/metadata.tsv", sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE) 
meta = metadata[,c("File.accession","Experiment.target","Assembly","Biological.replicate.s.")] %>%
  dplyr::filter(Assembly == "GRCh38")
colnames(meta) = c("file_name", "target_name", "assembly", "replicate")
meta = tbl_df(meta) %>% dplyr::mutate(sample_name = paste(target_name, replicate, sep = "_")) %>%
  tidyr::separate(target_name,c("factor_name", "species"), sep = "-")

#Import factor functions
factor_data = readr::read_csv("processed/annotations/fgwas/raw/RBP_summaries.csv")
factor_data = factor_data[,c("X1", "Splicing regulation", "3' end processing")]
colnames(factor_data) = c("factor_name","splicing", "three_end")

#Identify different classes of factors
three_end_only = dplyr::semi_join(meta, dplyr::filter(factor_data, three_end == 1, splicing == 0), by = "factor_name")
both = dplyr::semi_join(meta, dplyr::filter(factor_data, three_end == 1, splicing == 1), by = "factor_name")
splicing_only = dplyr::semi_join(meta, dplyr::filter(factor_data, three_end == 0, splicing == 1), by = "factor_name")

#Import peaks
peaks = loadNarrowPeaks("processed/annotations/fgwas/raw/eCLIP/", sample_names = meta$File.accession, 
                peaks_suffix = ".bed.gz", sub_dir = FALSE)
names(peaks) = meta$sample_name

#Three end peaks
peaks_twice_found = filterOverlaps(peaks[three_end_only$sample_name], minOverlapCount = 2)
three_end_peaks = listUnion(peaks_twice_found)
three_end_peaks_clean = keepSeqlevels(three_end_peaks, valid_chrs, pruning.mode = "tidy") %>% sort()
saveRDS(three_end_peaks_clean, "processed/annotations/fgwas/eCLIP_three_end_factors.rds")

#Three end peaks
peaks_twice_found = filterOverlaps(peaks[splicing_only$sample_name], minOverlapCount = 2)
splicing_peaks = listUnion(peaks_twice_found)
splicing_peaks_clean = keepSeqlevels(splicing_peaks, valid_chrs, pruning.mode = "tidy") %>% sort()
saveRDS(splicing_peaks_clean, "processed/annotations/fgwas/eCLIP_splicing_factors.rds")

peaks_twice_found = filterOverlaps(peaks[both$sample_name], minOverlapCount = 2)
both_peaks = listUnion(peaks_twice_found)
both_peaks_clean = keepSeqlevels(both_peaks, valid_chrs, pruning.mode = "tidy") %>% sort()
saveRDS(both_peaks_clean, "processed/annotations/fgwas/eCLIP_both_factors.rds")

