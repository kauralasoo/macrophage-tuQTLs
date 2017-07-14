library("dplyr")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Functions
calculatePositionPi1 <- function(p_matrix, p_threshold = 1e-5){
  
  up_down = dplyr::filter(p_matrix, upstream < p_threshold, !is.na(downstream))
  p1 = 1 - qvalue::qvalue(up_down$downstream)$pi0
  
  up_contained = dplyr::filter(p_matrix, upstream < p_threshold, !is.na(contained))
  p2 = 1 - qvalue::qvalue(up_contained$contained)$pi0
  
  down_contained = dplyr::filter(p_matrix, downstream < p_threshold, !is.na(contained))
  p3 = 1 - qvalue::qvalue(down_contained$contained)$pi0
  
  down_up = dplyr::filter(p_matrix, downstream < p_threshold, !is.na(upstream))
  p4 = 1 - qvalue::qvalue(down_up$upstream)$pi0
  
  contained_up = dplyr::filter(p_matrix, contained < p_threshold, !is.na(upstream))
  p5 = 1 - qvalue::qvalue(contained_up$upstream)$pi0
  
  contained_down = dplyr::filter(p_matrix, contained < p_threshold, !is.na(downstream))
  p6 = 1 - qvalue::qvalue(contained_down$downstream)$pi0
  
  #Combine into df
  result = dplyr::data_frame(first = c("start","start","end","end","mid","mid"),
                             second = c("end","mid","mid","start","start","end"),
                             pi1 = c(p1,p2,p3,p4,p5,p6))
  return(result)
}

calculateLeafCutterPi1 <- function(pvalues_df, leafcutter_names, p_threshold, max_cluster_count = 4){
  
  #Identify genes with multiple leafcutter clusters
  qtls = pvalues_df %>%
    dplyr::left_join(leafcutter_names, by = "phenotype_id") %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::arrange(gene_id, phenotype_id) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(event_count = length(gene_id)) %>%
    dplyr::filter(event_count > 1) %>%
    dplyr::filter(event_count <= max_cluster_count) %>%
    dplyr::arrange(gene_id) %>%
    dplyr::mutate(event_index = row_number()) %>%
    dplyr::ungroup()
  
  #Construct pairs of intron groups within genes
  group_pairs = dplyr::select(qtls, gene_id, group_id) %>% 
    dplyr::group_by(gene_id) %>% 
    purrrlyr::by_slice(~constructCombinatorialPairs(.$group_id), .collate = "rows")
  
  #Add first and second p-values
  p_df = dplyr::left_join(group_pairs, dplyr::transmute(qtls, first = group_id, p_beta_first = p_beta), by = "first") %>%
    dplyr::left_join(dplyr::transmute(qtls, second = group_id, p_beta_second = p_beta), by = "second")
  
  lc_pi1_first = 1-qvalue::qvalue(dplyr::filter(p_df, p_beta_first < p_threshold)$p_beta_second)$pi0
  lc_pi1_second = 1-qvalue::qvalue(dplyr::filter(p_df, p_beta_second < p_threshold)$p_beta_first)$pi0
  pi1_df = dplyr::data_frame(comparison = c("first_second", "second_first"), pi1 = c(lc_pi1_first, lc_pi1_second))
  
  return(pi1_df)
}

calculateConditionsPi1 <- function(pvalues_list, p_threshold = 1e-5){
  
  #Find all hits
  all_hits = purrr::map(pvalues_list, ~dplyr::filter(.,p_beta < p_threshold))
  
  #Calculate all pi1
  others = purrr::map(pvalues_list, ~dplyr::semi_join(., all_hits$naive, by = "group_id"))
  naive_pi1 = purrr::map_df(others[2:4], ~data_frame(first = "naive", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")
  
  others = purrr::map(pvalues_list, ~dplyr::semi_join(., all_hits$IFNg, by = "group_id"))
  ifng_pi1 = purrr::map_df(others[c(1,3,4)], ~data_frame(first = "IFNg", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")
  
  others = purrr::map(pvalues_list, ~dplyr::semi_join(., all_hits$SL1344, by = "group_id"))
  sl1344_pi1 = purrr::map_df(others[c(1,2,4)], ~data_frame(first = "SL1344", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")
  
  others = purrr::map(pvalues_list, ~dplyr::semi_join(., all_hits$IFNg_SL1344, by = "group_id"))
  IFNg_SL1344_pi1 = purrr::map_df(others[c(1,2,3)], ~data_frame(first = "IFNg_SL1344", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")
  
  condition_pi1_df = dplyr::bind_rows(naive_pi1, ifng_pi1, sl1344_pi1, IFNg_SL1344_pi1)
  return(condition_pi1_df)
}

constructCombinatorialPairs <- function(name_vector){
  pairs = combn(name_vector,2) %>% 
    t() %>% 
    dplyr::tbl_df()
  colnames(pairs) = c("first", "second")
  return(pairs)
}

#Import SummarizedExperiments
se_revised = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds")
se_leafcutter = readRDS("results/SummarizedExperiments/salmonella_leafcutter_counts.rds")
se_ensembl = readRDS("results/SummarizedExperiments/salmonella_salmon_Ensembl_87.rds")

#Construct gene name map
revised_names = tbl_df2(rowData(se_revised)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id = ensembl_gene_id)
leafcutter_names = tbl_df2(rowData(se_leafcutter)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id = ensembl_gene_id)
ensembl_names = tbl_df2(rowData(se_ensembl)) %>%
  dplyr::transmute(phenotype_id = transcript_id, gene_name, gene_id)

#Import QTLs
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")


#### reviseAnnotations ####
revised_qtls = salmonella_qtls$reviseAnnotations_groupwise[1:4]

#Add names to QTLs
named_qtls = purrr::map(revised_qtls, ~dplyr::left_join(.,revised_names, by = "phenotype_id") %>%
                        tidyr::separate("group_id", c("gene1","group_name","position"), sep = "\\.", remove = FALSE) %>%
                        dplyr::select(-gene1))

#Identify best group in each condition
#Count QTLs per group
best_group = purrr::map(named_qtls, ~dplyr::group_by(., gene_id, gene_name, group_name) %>% 
  dplyr::arrange(phenotype_id) %>%
  dplyr::filter(p_fdr < 0.1) %>% 
  dplyr::summarize(qtl_count = length(phenotype_id), min_p = min(p_beta)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::arrange(gene_id, -qtl_count) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::ungroup())

#Filter best groups in each condition
best_filtered = purrr::map2(named_qtls, best_group, ~dplyr::semi_join(.x,.y, by = c("gene_id","group_name")))

#### Pi1 analysis between positions ####
#Contruct a matrix of p-values for each position
p_matrix_list = purrr::map(best_filtered, ~dplyr::select(., gene_id, gene_name, position, p_beta) %>% 
  tidyr::spread(position, p_beta))

#Calculate pi1 between all positions in each condition
position_pi1_df = purrr::map_df(p_matrix_list, ~calculatePositionPi1(.,p_threshold = 0.001), .id = "condition_name") %>%
  dplyr::mutate(comparison = paste(first, second, sep = "_"))

#Make a boxplot
pi1_positions_plot = ggplot(position_pi1_df, aes(x = comparison, y = pi1)) + 
  geom_boxplot() +
  geom_point() + 
  scale_y_continuous(limits = c(0,1))
ggsave("results/figures/pi1_position_plot.pdf",plot = pi1_positions_plot, width = 5, height = 4)


#### leafcutter ####
leafcutter_pi1_df = purrr::map_df(salmonella_qtls$leafcutter[1:4], 
                                  ~calculateLeafCutterPi1(., leafcutter_names, p_threshold = 0.001, max_cluster_count = 4), .id = "condition_name")
pi1_leafcutter_plot = ggplot(leafcutter_pi1_df, aes(x = "Pi1", y = pi1)) + 
  geom_boxplot() +
  geom_point() + 
  scale_y_continuous(limits = c(0,1))
ggsave("results/figures/pi1_leafcutter_plot.pdf",plot = pi1_leafcutter_plot, width = 2, height = 4)




#### Pi1 analysis between conditions ####

# Calculate all pi1 values
fc = calculateConditionsPi1(salmonella_qtls$featureCounts[1:4], p_threshold = 0.001)
re = calculateConditionsPi1(salmonella_qtls$reviseAnnotations[1:4], p_threshold = 0.001)
lc = calculateConditionsPi1(salmonella_qtls$leafcutter[1:4], p_threshold = 0.001)
en = calculateConditionsPi1(salmonella_qtls$Ensembl_87[1:4], p_threshold = 0.001)
pi_df = purrr::map_df(list(featureCounts = fc, reviseAnnotations = re, leafcutter = lc, Ensembl_87 = en), identity, .id = "method")


#Make a boxplot
pi1_condition_plot = ggplot(pi_df, aes(x = method, y = pi1)) + 
  geom_boxplot() +
  geom_point() + 
  scale_y_continuous(limits = c(0,1))
ggsave("results/figures/pi1_condition_plot.pdf",plot = pi1_condition_plot, width = 5, height = 4)





