library("readr")
library("tidyr")
library("dplyr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")


#Import caQTL lead variants
caQTLs = readr::read_tsv("processed/annotations/reference_QTLs/caQTL_naive_RASQUAL_lead_only.txt.gz") %>%
  dplyr::filter(p_eigen < fdr_thresh)

#Import eQTL summaries
tuQTLs = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")$reviseAnnotations$naive %>%
  dplyr::filter(p_fdr < 0.1) %>%
  tidyr::separate(group_id, c("gene_id", "position"), sep = "\\.", remove = FALSE)

#Import txrevise metadata
revised_event_meta = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2()

#Import genotypes
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")

#Identify trait-SNP pairs
tuQTL_pairs = dplyr::transmute(tuQTLs, gene_id = phenotype_id, snp_id) %>%
  addVariantCoords(., vcf_file$snpspos)
caQTL_pairs = dplyr::transmute(caQTLs, peak_id = gene_id, snp_id) %>%
  addVariantCoords(., vcf_file$snpspos)
  
#Find overlaps chr by chr
chr_list = idVectorToList(unique(tuQTL_pairs$chr))
overlap_list = purrr::map(chr_list, ~findGWASOverlaps(dplyr::filter(tuQTL_pairs, chr == .) %>% 
                                                        dplyr::select(gene_id, snp_id), 
                                                      dplyr::filter(caQTL_pairs, chr == .), 
                                                      vcf_file, max_distance = 5e5, min_r2 = 0.8))
tuQTL_caQTL_overlaps = purrr::map_df(overlap_list, identity) %>% 
  dplyr::rename(phenotype_id = gene_id) %>%
  dplyr::left_join(dplyr::select(tuQTLs, phenotype_id, position))

overlap = dplyr::group_by(tuQTL_caQTL_overlaps, position) %>% 
  dplyr::select(phenotype_id) %>%
  dplyr::distinct() %>%
  dplyr::summarise(overlap_count = length(phenotype_id))

all = dplyr::group_by(tuQTLs, position) %>% dplyr::select(phenotype_id) %>%
  dplyr::distinct() %>%
  dplyr::summarise(total_count = length(phenotype_id))

fisher.test(matrix(c(194, 993-194, 248, 1398-248), ncol = 2, byrow = T))
fisher.test(matrix(c(194, 993-194, 176, 1105-176), ncol = 2, byrow = T))

fisher.test(matrix(c(62, 221-62, 248, 1398-248), ncol = 2, byrow = T))


#restrict analysis to true promoter events
promoter_events = readRDS("results/reviseAnnotations/promoter_events.rds")
true_promoters = dplyr::filter(promoter_events, contained == 0) %>%
  dplyr::transmute(group_id, is_promoter = TRUE) %>%
  dplyr::left_join(dplyr::transmute(revised_event_meta, group_id, phenotype_id = transcript_id))
  
dplyr::semi_join(tuQTLs, true_promoters, by = "phenotype_id")
dplyr::semi_join(tuQTL_caQTL_overlaps, true_promoters, by = "phenotype_id") %>%
  dplyr::select(phenotype_id) %>%
  dplyr::distinct()




