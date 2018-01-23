library("data.table")

qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

#Ends
new_end_qtls = purrr::map_df(qtls$txrevise_ends, ~dplyr::filter(.,p_fdr < 0.1) %>% nrow())
orig_end_qtls = purrr::map_df(qtls$reviseAnnotations, ~dplyr::filter(.,group_id %like% "downstream") %>% 
                dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>% 
                dplyr::filter(p_fdr < 0.1) %>% nrow())

new_end_qtls/orig_end_qtls


#Ends
new_promoters_qtls = purrr::map_df(qtls$txrevise_promoters, ~dplyr::filter(.,p_fdr < 0.1) %>% nrow())
orig_promoters_qtls = purrr::map_df(qtls$reviseAnnotations, ~dplyr::filter(.,group_id %like% "upstream") %>% 
                                dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>% 
                                dplyr::filter(p_fdr < 0.1) %>% nrow())
new_promoters_qtls/orig_promoters_qtls