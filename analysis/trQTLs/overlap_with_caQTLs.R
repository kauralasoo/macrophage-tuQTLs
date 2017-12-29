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
                                                      vcf_file, max_distance = 5e5, min_r2 = 0.9))
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

all_overlaps = dplyr::left_join(overlap, all) %>% 
  dplyr::mutate(fraction = overlap_count/total_count)

#restrict analysis to true promoter events
puQTLs = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")$txrevise_promoters$naive %>%
  dplyr::filter(p_fdr < 0.1) %>%
  tidyr::separate(group_id, c("gene_id", "position"), sep = "\\.", remove = FALSE)

#Identify trait-SNP pairs
puQTL_pairs = dplyr::transmute(puQTLs, gene_id = phenotype_id, snp_id) %>%
  addVariantCoords(., vcf_file$snpspos)

#Find overlaps chr by chr
chr_list = idVectorToList(unique(puQTL_pairs$chr))
overlap_list = purrr::map(chr_list, ~findGWASOverlaps(dplyr::filter(puQTL_pairs, chr == .) %>% 
                                                        dplyr::select(gene_id, snp_id), 
                                                      dplyr::filter(caQTL_pairs, chr == .), 
                                                      vcf_file, max_distance = 5e5, min_r2 = 0.9))
tuQTL_caQTL_overlaps = purrr::map_df(overlap_list, identity) %>% 
  dplyr::rename(phenotype_id = gene_id) %>%
  dplyr::left_join(dplyr::select(puQTLs, phenotype_id, position))

overlap = dplyr::group_by(tuQTL_caQTL_overlaps, position) %>% 
  dplyr::select(phenotype_id) %>%
  dplyr::distinct() %>%
  dplyr::summarise(overlap_count = length(phenotype_id))

all = dplyr::group_by(puQTLs, position) %>% dplyr::select(phenotype_id) %>%
  dplyr::distinct() %>%
  dplyr::summarise(total_count = length(phenotype_id))

promoter_overlaps = dplyr::left_join(overlap, all) %>% 
  dplyr::mutate(fraction = overlap_count/total_count)

#Save all ratios to disk
caQTL_olaps = dplyr::bind_rows(all_overlaps, promoter_overlaps)
write.table(caQTL_olaps, "results/tables/caQTL_overlaps.txt", sep = "\t", quote = F, row.names = F)

#R2 > 0.8
fisher.test(matrix(c(167, 786-167, 248, 1398-248), ncol = 2, byrow = T))
fisher.test(matrix(c(167, 786-167, 176, 1105-176), ncol = 2, byrow = T))
fisher.test(matrix(c(167, 786-167, (248+176), (1398+1105)-(248+176)), ncol = 2, byrow = T))

#R2 > 0.9
fisher.test(matrix(c(124, 786-124, 152, 1398-152), ncol = 2, byrow = T))
fisher.test(matrix(c(124, 786-124, 103, 1105-103), ncol = 2, byrow = T))
fisher.test(matrix(c(124, 786-124, (152+103), (1398+1105)-(152+103)), ncol = 2, byrow = T))




