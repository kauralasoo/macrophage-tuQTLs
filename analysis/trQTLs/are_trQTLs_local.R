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
position_pi1_df = purrr::map_df(p_matrix_list, calculatePositionPi1, .id = "condition_name") %>%
  dplyr::mutate(comparison = paste(first, second, sep = "_"))

#Make a boxplot
pi1_positions_plot = ggplot(position_pi1_df, aes(x = comparison, y = pi1)) + 
  geom_boxplot() +
  geom_point() + 
  scale_y_continuous(limits = c(0,1))
ggsave("results/figures/pi1_position_plot.pdf",plot = pi1_positions_plot, width = 5, height = 4)


#### Pi1 analysis between conditions ####
all_hits = purrr::map(best_filtered, ~dplyr::filter(.,p_beta < 1e-5))
others = purrr::map(named_qtls, ~dplyr::semi_join(., all_hits$naive, by = "group_id"))
naive_pi1 = purrr::map_df(others[2:4], ~data_frame(first = "naive", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")

others = purrr::map(named_qtls, ~dplyr::semi_join(., all_hits$IFNg, by = "group_id"))
ifng_pi1 = purrr::map_df(others[c(1,3,4)], ~data_frame(first = "IFNg", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")

others = purrr::map(named_qtls, ~dplyr::semi_join(., all_hits$SL1344, by = "group_id"))
sl1344_pi1 = purrr::map_df(others[c(1,2,4)], ~data_frame(first = "SL1344", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")

others = purrr::map(named_qtls, ~dplyr::semi_join(., all_hits$IFNg_SL1344, by = "group_id"))
IFNg_SL1344_pi1 = purrr::map_df(others[c(1,2,3)], ~data_frame(first = "IFNg_SL1344", pi1 = 1 - qvalue::qvalue(.$p_beta)$pi0), .id = "second")

condition_pi1_df = dplyr::bind_rows(naive_pi1, ifng_pi1, sl1344_pi1, IFNg_SL1344_pi1)




