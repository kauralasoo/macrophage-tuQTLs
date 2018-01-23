library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import gene metadata
revised_gene_metadata = readRDS("results/SummarizedExperiments/salmonella_salmon_reviseAnnotations.rds") %>%
  rowData(.) %>% tbl_df2()
group_map = dplyr::select(revised_gene_metadata, transcript_id, group_id) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(transcript_count = length(transcript_id)) %>%
  dplyr::ungroup()

#Import QTLs
salmonella_df = readRDS("results/trQTLs/variance_explained/salmonella_compiled_varExp.rds")
acldl_df = readRDS("results/trQTLs/variance_explained/acldl_compiled_varExp.rds")
qtls_df = dplyr::bind_rows(salmonella_df, acldl_df) %>%
  dplyr::left_join(group_map, by = c("phenotype_id" = "transcript_id")) %>%
  #dplyr::left_join(true_promoters, by = "group_id") %>%
  #dplyr::mutate(is_promoter = ifelse(is.na(is_promoter), FALSE, is_promoter)) %>%
  dplyr::mutate(is_response = ifelse(interaction_fraction > 0.5 & p_fdr < 0.1, TRUE, FALSE))

txrevise_events = dplyr::filter(qtls_df, quant == "reviseAnnotations") %>%
  tidyr::separate(phenotype_id, c("ensembl_gene_id","grp", "position", "transcript_id"), sep = "\\.", remove = FALSE)

promoter_events = dplyr::filter(qtls_df, quant == "txrevise_promoters") %>%
  tidyr::separate(phenotype_id, c("ensembl_gene_id","grp", "position", "transcript_id"), sep = "\\.", remove = FALSE)

end_events = dplyr::filter(qtls_df, quant == "txrevise_ends") %>%
  tidyr::separate(phenotype_id, c("ensembl_gene_id","grp", "position", "transcript_id"), sep = "\\.", remove = FALSE)

#Count all qtls
other_counts = txrevise_events %>% 
  group_by(condition, position) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), response_count = sum(is_response, na.rm = T)) %>% 
  dplyr::mutate(response_fraction = response_count/qtl_count)

#Event names
event_names = data_frame(position = c("upstream", "contained","downstream"), 
                         event_type = factor(c("promoters", "middle exons", "3' ends"), levels = c("promoters", "middle exons", "3' ends")))

count_df = other_counts %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(event_names)

promoter_plot = ggplot(count_df, aes(x = event_type, y = response_fraction, group = figure_name, color = figure_name)) + 
  geom_point() + 
  geom_line() + 
  theme_light() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
  ylab("Response QTL fraction") +
  scale_color_manual(name = "condition", values = conditionPalette()) + 
  coord_cartesian(ylim = c(0,0.17))
ggsave("results/figures/response_fraction_by_event_type.pdf", plot = promoter_plot, height = 2.5, width = 3)


#### Look at true promoter events ####
promoter_counts = promoter_events %>% 
  group_by(condition, position) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), response_count = sum(is_response, na.rm = T)) %>% 
  dplyr::mutate(response_fraction = response_count/qtl_count)

end_counts = end_events %>% 
  group_by(condition, position) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), response_count = sum(is_response, na.rm = T)) %>% 
  dplyr::mutate(response_fraction = response_count/qtl_count)

combined_counts = dplyr::bind_rows(dplyr::filter(other_counts, position == "contained"), promoter_counts, end_counts)

event_names = data_frame(position = c("upstream", "contained","downstream"), event_type = factor(c("promoters", "middle exons", "3' ends"), levels = c("promoters", "middle exons", "3' ends")))

count_df = combined_counts %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(event_names)
write.table(count_df, "results/tables/promoter_response_QTL_fraction.txt", sep = "\t", quote = F, row.names = F)

promoter_plot = ggplot(count_df, aes(x = event_type, y = response_fraction, group = figure_name, color = figure_name)) + 
  geom_point() + 
  geom_line() + 
  theme_light() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
  ylab("Response QTL fraction") +
  scale_color_manual(name = "condition", values = conditionPalette()) +
  coord_cartesian(ylim = c(0,0.17))
ggsave("results/figures/response_fraction_by_promoter.pdf", plot = promoter_plot, height = 2.5, width = 3)

#Perform Fisher's tests
fisher.test(matrix(c(84, 490-84, 104, 843-104), ncol = 2, byrow = T))
fisher.test(matrix(c(79, 743-79, 75, 1273-75), ncol = 2, byrow = T))
fisher.test(matrix(c(72, 580-72, 106, 1049-106), ncol = 2, byrow = T))

#Combined:
#IFNg + Salmonella
fisher.test(matrix(c(84, 490-84, 104+75, (843+664)-(104+75)), ncol = 2, byrow = T))

#IFNg
fisher.test(matrix(c(79, 743-79, 60+75, (1273+1025)-(60+75)), ncol = 2, byrow = T))

#Salmonella
fisher.test(matrix(c(72, 580-72, 106+75, (1049+851)-(106+75)), ncol = 2, byrow = T))

#AcLDL
fisher.test(matrix(c(24, 586-24, 45+19, (1076+802)-(45+19)), ncol = 2, byrow = T))

#Combined p-value
metap::sumlog(c(0.4446, 0.04989, 2.302e-05, 0.003439))

metap::sumlog(c(0.4446, 0.04989, 2.302e-05, 0.003439))
