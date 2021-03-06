
#' Add ATAC peak annotations to the eQTL credible set
#'
#' @param eqtl_credible_set 
#' @param atac_peak_metadata 
#' @param atac_tabix_path 
#'
#' @return Credible set with ATAC peak annotations
#' @export
annotateCredibleSet <- function(eqtl_credible_set, atac_peak_metadata, atac_tabix_path){
  
  #Construct Granges for the credible set SNPs
  variant_ranges = dplyr::transmute(eqtl_credible_set, seqnames = chr, start = pos, end = pos, strand = "+", snp_id) %>%
    dataFrameToGRanges()
  
  #Add overlapping peaks
  peak_finemapping = addOverlappingPeaks(eqtl_credible_set, atac_peak_metadata)
  
  #Find most associated peaks for each SNP
  associated_peaks = tabixFetchSNPs(variant_ranges, atac_tabix_path)
  if (!is.null(associated_peaks)){
    associated_peaks = associated_peaks %>% 
      dplyr::group_by(snp_id) %>% 
      dplyr::arrange(p_nominal) %>% 
      dplyr::filter(row_number() == 1) %>% 
      dplyr::ungroup() %>% 
      dplyr::transmute(assoc_peak_id = gene_id, snp_id, p_nominal_peak = p_nominal)
    
    #Add associated peaks
    peak_finemapping = dplyr::left_join(peak_finemapping, associated_peaks, by = "snp_id")
  }else{
    peak_finemapping = dplyr::mutate(peak_finemapping, assoc_peak_id = NA)
  }
  
  return(peak_finemapping)
}


#Calculate summary stats on the credible sets
summariseCredibleSets <- function(credible_sets_df){
  
  #Calculate credible set size
  credible_set_size = dplyr::group_by(credible_sets_df, gene_id) %>% 
    dplyr::arrange(p_nominal) %>% 
    dplyr::summarise(credible_set_size = length(snp_id))
  
  #Remove SNPs that do not overlap peaks
  peaks_filtered = dplyr::filter(credible_sets_df, !is.na(overlap_peak_id))
  
  #Find peaks where at least one variant lies within the peak
  overlap_snp_count = dplyr::filter(peaks_filtered, gene_id == overlap_peak_id) %>% 
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(overlap_snp_count = length(snp_id)) %>%
    dplyr::ungroup()
  
  #Count the number of peaks overlapping SNPs in the credible set
  overlap_peak_count = peaks_filtered %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::summarise(overlap_peak_count = length(overlap_peak_id)) %>% 
    ungroup()
  
  #Find lead SNPs
  lead_snps = dplyr::group_by(peak_cs_df, gene_id) %>% dplyr::arrange(p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>% dplyr::ungroup() %>%
    dplyr::left_join(credible_set_size, by = "gene_id") %>%
    dplyr::left_join(overlap_snp_count, by = "gene_id") %>%
    dplyr::mutate(overlap_snp_count = ifelse(is.na(overlap_snp_count), 0, overlap_snp_count)) %>%
    dplyr::left_join(overlap_peak_count, by = "gene_id") %>%
    dplyr::mutate(overlap_peak_count = ifelse(is.na(overlap_peak_count), 0, overlap_peak_count))
  
  return(lead_snps)
}

constructClustersFromGenePairs <- function(gene_pairs, cluster_name_prefix = "gene_cluster_"){
  assertthat::assert_that(ncol(gene_pairs) == 2)
  
  #Construct igraph object and extract clusters
  graph = igraph::graph_from_data_frame(gene_pairs, directed = FALSE)
  clusters = igraph::clusters(graph)
  
  #Make df
  clusters_df = data_frame(gene_id = names(clusters$membership), cluster_number = clusters$membership) %>% 
    arrange(cluster_number) %>% 
    dplyr::mutate(cluster_id = paste0(cluster_name_prefix, cluster_number)) %>% 
    dplyr::select(gene_id, cluster_id)
  
  #Count genes in clusters
  gene_counts = dplyr::group_by(clusters_df, cluster_id) %>% 
    dplyr::summarise(gene_count = length(gene_id))
  clusters_df = dplyr::left_join(clusters_df, gene_counts, by = "cluster_id")
  return(clusters_df)
  
}

credibleSetJaccard <- function(row_df, master_credible_sets, dependent_credible_sets){
  
  #Extract credible sets
  master_cs = master_credible_sets[[row_df$master_condition]][[row_df$master_id]] %>% dplyr::arrange(p_nominal)
  dependent_cs = dependent_credible_sets[[row_df$dependent_condition]][[row_df$dependent_id]] %>% dplyr::arrange(p_nominal)

  #Calculate two types of jaccard indexes
  jaccard = length(intersect(master_cs$snp_id, dependent_cs$snp_id)) / length(union(master_cs$snp_id, dependent_cs$snp_id))
  min_jaccard = length(intersect(master_cs$snp_id, dependent_cs$snp_id)) / 
    min(length(master_cs$snp_id), length(dependent_cs$snp_id))
  res = data_frame(jaccard = jaccard, min_jaccard = min_jaccard, snp_id = master_cs$snp_id[1], p_nominal = master_cs$p_nominal[1])
  return(res)
}

#Convert credible sets into gigantic data frame
credibleSetsToDf <- function(credible_sets){
  result = purrr::map_df(credible_sets, ~purrr::map_df(., 
                  ~dplyr::mutate(.,chr = as.character(chr))) %>% 
                  dplyr::filter(chr != "X"), .id = "condition_name")
  return(result)
}

credibleSetsToGranges <- function(credible_sets_df){
  granges = credible_sets_df %>% 
    dplyr::transmute(condition_name, gene_id, snp_id, chr, seqnames = chr, start = pos, end = pos, strand = "+") %>% 
    dataFrameToGRanges()
  return(granges)
}

findCredibleSetOverlaps <- function(master_credible_set, dependent_credible_set){
  #Find overlaps
  olaps = findOverlaps(master_credible_set, dependent_credible_set)
  
  #Eextract query overlaps
  queries = master_credible_set[queryHits(olaps),] %>% 
    elementMetadata() %>% as.data.frame() %>% 
    dplyr::transmute(master_condition = condition_name, master_id = gene_id, chr1 = chr) %>% 
    tbl_df()
  #Extract subject overlaps
  subjects = dependent_credible_set[subjectHits(olaps),] %>% 
    elementMetadata() %>% as.data.frame() %>% 
    dplyr::transmute(dependent_condition = condition_name, dependent_id = gene_id, chr2 = chr) %>% 
    tbl_df()
  
  #Identify all genes that share at least one SNP in their credible set
  pairwise_shared = dplyr::bind_cols(queries, subjects) %>% 
    unique() %>% 
    dplyr::filter(chr1 == chr2) %>% 
    dplyr::filter(master_id != dependent_id) %>% 
    dplyr::select(master_id, master_condition, dependent_id, dependent_condition) %>% 
    unique()
  return(pairwise_shared)
}
