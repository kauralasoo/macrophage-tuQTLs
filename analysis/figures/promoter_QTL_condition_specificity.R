library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

#Import true promoter events
promoter_events = readRDS("results/reviseAnnotations/promoter_events.rds")
true_promoters = dplyr::filter(promoter_events, contained == 0) %>%
  dplyr::transmute(group_id, is_promoter = TRUE)

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
  dplyr::left_join(true_promoters, by = "group_id") %>%
  dplyr::mutate(is_promoter = ifelse(is.na(is_promoter), FALSE, is_promoter)) %>%
  dplyr::mutate(is_response = ifelse(interaction_fraction > 0.5 & p_fdr < 0.1, TRUE, FALSE)) %>%
  dplyr::filter(quant == "reviseAnnotations") %>%
  dplyr::filter(transcript_count <= 10) %>%
  tidyr::separate(phenotype_id, c("ensembl_gene_id","grp", "position", "transcript_id"), sep = "\\.", remove = FALSE)

#Count promoter qtls
other_counts = qtls_df %>% 
  group_by(condition, position) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), response_count = sum(is_response, na.rm = T)) %>% 
  dplyr::mutate(response_fraction = response_count/qtl_count)

#Event names
event_names = data_frame(position = c("upstream", "contained","downstream"), event_type = factor(c("start", "middle", "end"), levels = c("start", "middle", "end")))

count_df = other_counts %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(event_names)

promoter_plot = ggplot(count_df, aes(x = event_type, y = response_fraction, group = figure_name, color = figure_name)) + 
  geom_point() + 
  geom_line() + 
  theme_light() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
  ylab("Response QTL fraction")
ggsave("results/figures/response_fraction_by_event_type.pdf", plot = promoter_plot, height = 2.5, width = 3)


#### Look at true promoter events ####
promoter_counts = dplyr::filter(qtls_df, is_promoter) %>% 
  group_by(condition, position) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), response_count = sum(is_response, na.rm = T)) %>% 
  dplyr::mutate(response_fraction = response_count/qtl_count)
combined_counts = dplyr::bind_rows(dplyr::filter(other_counts, position != "upstream"), promoter_counts)

event_names = data_frame(position = c("upstream", "contained","downstream"), event_type = factor(c("promoter", "middle", "end"), levels = c("promoter", "middle", "end")))

count_df = combined_counts %>%
  dplyr::rename(condition_name = condition) %>%
  dplyr::left_join(conditionFriendlyNames()) %>%
  dplyr::left_join(event_names)

promoter_plot = ggplot(count_df, aes(x = event_type, y = response_fraction, group = figure_name, color = figure_name)) + 
  geom_point() + 
  geom_line() + 
  theme_light() +
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
  ylab("Response QTL fraction")
ggsave("results/figures/response_fraction_by_promoter.pdf", plot = promoter_plot, height = 2.5, width = 3)


promoter_counts = dplyr::filter(qtls_df, position == "upstream") %>% 
  group_by(condition, is_promoter) %>% 
  dplyr::summarise(qtl_count = length(phenotype_id), response_count = sum(is_response, na.rm = T)) %>% 
  dplyr::mutate(response_fraction = response_count/qtl_count)
