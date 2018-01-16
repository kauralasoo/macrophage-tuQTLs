library("dplyr")
library("readr")
library("rtracklayer")
library("tximport")
library("devtools")
library("GenomicFeatures")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")


convertChrNames <- function(granges, to = c("UCSC", "Ensembl")) {
  ensembl <- c("1", "2", "3", "4", "5", "6", "7", "X", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "20", "Y", "19", "22", "21", "GL000192.1", "GL000225.1", "GL000194.1", "GL000193.1", "GL000200.1", "GL000222.1", "GL000212.1", "GL000195.1", "GL000223.1", "GL000224.1", "GL000219.1", "GL000205.1", "GL000215.1", "GL000216.1", "GL000217.1", "GL000199.1", "GL000211.1", "GL000213.1", "GL000220.1", "GL000218.1", "GL000209.1", "GL000221.1", "GL000214.1", "GL000228.1", "GL000227.1", "GL000191.1", "GL000208.1", "GL000198.1", "GL000204.1", "GL000233.1", "GL000237.1", "GL000230.1", "GL000242.1", "GL000243.1", "GL000241.1", "GL000236.1", "GL000240.1", "GL000206.1", "GL000232.1", "GL000234.1", "GL000202.1", "GL000238.1", "GL000244.1", "GL000248.1", "GL000196.1", "GL000249.1", "GL000246.1", "GL000203.1", "GL000197.1", "GL000245.1", "GL000247.1", "GL000201.1", "GL000235.1", "GL000239.1", "GL000210.1", "GL000231.1", "GL000229.1", "MT", "GL000226.1", "GL000207.1", "HSCHR4_1", "HSCHR17_1", "HSCHR6_MHC_SSTO", "HSCHR6_MHC_QBL", "HSCHR6_MHC_MCF", "HSCHR6_MHC_MANN", "HSCHR6_MHC_COX", "HSCHR6_MHC_DBB", "HSCHR6_MHC_APD")
  ucsc <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chrY", "chr19", "chr22", "chr21", "chr1_gl000192_random", "chrUn_gl000225", "chr4_gl000194_random", "chr4_gl000193_random", "chr9_gl000200_random", "chrUn_gl000222", "chrUn_gl000212", "chr7_gl000195_random", "chrUn_gl000223", "chrUn_gl000224", "chrUn_gl000219", "chr17_gl000205_random", "chrUn_gl000215", "chrUn_gl000216", "chrUn_gl000217", "chr9_gl000199_random", "chrUn_gl000211", "chrUn_gl000213", "chrUn_gl000220", "chrUn_gl000218", "chr19_gl000209_random", "chrUn_gl000221", "chrUn_gl000214", "chrUn_gl000228", "chrUn_gl000227", "chr1_gl000191_random", "chr19_gl000208_random", "chr9_gl000198_random", "chr17_gl000204_random", "chrUn_gl000233", "chrUn_gl000237", "chrUn_gl000230", "chrUn_gl000242", "chrUn_gl000243", "chrUn_gl000241", "chrUn_gl000236", "chrUn_gl000240", "chr17_gl000206_random", "chrUn_gl000232", "chrUn_gl000234", "chr11_gl000202_random", "chrUn_gl000238", "chrUn_gl000244", "chrUn_gl000248", "chr8_gl000196_random", "chrUn_gl000249", "chrUn_gl000246", "chr17_gl000203_random", "chr8_gl000197_random", "chrUn_gl000245", "chrUn_gl000247", "chr9_gl000201_random", "chrUn_gl000235", "chrUn_gl000239", "chr21_gl000210_random", "chrUn_gl000231", "chrUn_gl000229", "chrM", "chrUn_gl000226", "chr18_gl000207_random", "chr4_ctg9_hap1", "chr17_ctg5_hap1", "chr6_ssto_hap7", "chr6_qbl_hap6", "chr6_mcf_hap5", "chr6_mann_hap4", "chr6_cox_hap2", "chr6_dbb_hap3", "chr6_apd_hap1")
  if (to == "Ensembl") {
    seqlevels(granges) <- ensembl[match(seqlevels(granges), ucsc)]
  }
  else if (to == "UCSC") {
    seqlevels(granges) <- ucsc[match(seqlevels(granges), ucsc)]
  }
  granges
}
valid_chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chrY", "chr19", "chr22", "chr21")

#Import variant information
variant_information = importVariantInformation("results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz")
variant_ranges = GRanges(seqnames = variant_information$chr, ranges = IRanges(start = variant_information$pos, end = variant_information$pos), strand = "*")
elementMetadata(variant_ranges) = as.data.frame(dplyr::select(variant_information, snp_id, MAF, type))

#Import ATAC peaks
atac_peaks = rtracklayer::import.gff3("processed/annotations/reference_QTLs/ATAC_consensus_peaks.gff3")
olaps = findOverlaps(variant_ranges, atac_peaks, ignore.strand=TRUE)
variant_ranges$atac_peak = 0
variant_ranges[queryHits(olaps)]$atac_peak = 1

#Imort Poly(A) peaks
extra_cols = c(motif = "character", gene_name = "character")
polyA = rtracklayer::import.bed("processed/annotations/fgwas/raw/PolyA.hg38.bed.gz", extraCols = extra_cols)
polyA_clean = keepSeqlevels(polyA, valid_chrs, pruning.mode = "tidy") %>% 
  convertChrNames(to = "Ensembl")
polyA_wide_df = tbl_df(as.data.frame(polyA_clean)) %>% 
  dplyr::mutate(centre = floor(start+((end-start)/2))) %>%
  dplyr::mutate(new_start = centre-24, new_end = centre+25)
polyA_wide = GRanges(seqnames = polyA_wide_df$seqnames, 
                     IRanges(start = polyA_wide_df$new_start, end = polyA_wide_df$new_end), 
                     strand = polyA_wide_df$strand) %>% reduce()

olaps = findOverlaps(variant_ranges, polyA_wide, ignore.strand=TRUE)
variant_ranges$polyA_site = 0
variant_ranges[queryHits(olaps)]$polyA_site = 1

#Import gene annotations
txdb = loadDb("../../annotations/GRCh38/genes/Ensembl_87/TranscriptDb_GRCh38_87.db")
tx_meta = readRDS("../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.compiled_tx_metadata.rds")
protein_transcripts = dplyr::filter(tx_meta, gene_biotype == "protein_coding", transcript_biotype == "protein_coding")

#5'UTRs
fiveUTRs = fiveUTRsByTranscript(txdb, use.names = TRUE)
fiveUTRs = fiveUTRs[intersect(names(fiveUTRs), protein_transcripts$ensembl_transcript_id)]
olaps = findOverlaps(variant_ranges, fiveUTRs, ignore.strand=TRUE)
variant_ranges$fiveUTR = 0
variant_ranges[queryHits(olaps)]$fiveUTR = 1

#3'UTRs
threeUTRs = threeUTRsByTranscript(txdb, use.names = TRUE)
threeUTRs = threeUTRs[intersect(names(threeUTRs), protein_transcripts$ensembl_transcript_id)]
olaps = findOverlaps(variant_ranges, threeUTRs, ignore.strand=TRUE)
variant_ranges$threeUTR = 0
variant_ranges[queryHits(olaps)]$threeUTR = 1

#CDSs
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdss[intersect(names(cdss), protein_transcripts$ensembl_transcript_id)]
olaps = findOverlaps(variant_ranges, cdss, ignore.strand=TRUE)
variant_ranges$CDS = 0
variant_ranges[queryHits(olaps)]$CDS = 1

#introns
introns = intronsByTranscript(txdb, use.names = TRUE)
introns = cdss[intersect(names(introns), protein_transcripts$ensembl_transcript_id)]
olaps = findOverlaps(variant_ranges, introns, ignore.strand=TRUE)
variant_ranges$intron = 0
variant_ranges[queryHits(olaps)]$intron = 1

#Promoters
promoters = promoters(txdb, upstream = 2000, downstream = 200)
promoters = promoters[promoters$tx_name %in% protein_transcripts$ensembl_transcript_id,] %>% reduce()
olaps = findOverlaps(variant_ranges, promoters, ignore.strand=TRUE)
variant_ranges$promoter = 0
variant_ranges[queryHits(olaps)]$promoter = 1

#CTCF
ctcf = loadNarrowPeaks("processed/annotations/fgwas/raw/", sample_names = c("CTCF_MAC"), sub_dir = FALSE)[[1]]
olaps = findOverlaps(variant_ranges, ctcf, ignore.strand=TRUE)
variant_ranges$CTCF = 0
variant_ranges[queryHits(olaps)]$CTCF = 1

#Export fgwas annotations
df = as.data.frame(variant_ranges) %>% 
  dplyr::select(-width, -strand, -end) %>%
  dplyr::rename(chr = seqnames, pos = start)

gz1 <- gzfile("processed/annotations/fgwas/fgwas_annotations.txt.gz", "w")
write.table(df, gz1, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz1)
