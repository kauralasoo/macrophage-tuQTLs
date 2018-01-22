library("dplyr")
library("readr")
library("tximport")
library("devtools")
library("ggplot2")
library("SummarizedExperiment")
load_all("../seqUtils/")
load_all("analysis/housekeeping/")

file_list = list.files("processed/salmonella/fgwas/output/", recursive = T)

df = data_frame(file_path = file_list[grep("params", file_list)]) %>%
  tidyr::separate(file_path, c("quant", "cond_string"), sep ="\\/", remove = F) %>% 
  tidyr::separate(cond_string, c("condition_name","empty"), sep = ".params" ) %>%
  dplyr::select(-empty) %>%
  dplyr::mutate(full_path = file.path("processed/salmonella/fgwas/output", file_path))

# Import files
file_paths = setNames(df$full_path, df$file_path)
all_data = purrr::map_df(file_paths, ~readr::read_delim(., delim = " ", col_types = "cddd"), .id = "file_path") %>%
  dplyr::mutate(CI_hi = ifelse(is.na(CI_hi), estimate + (estimate-CI_lo), CI_hi)) %>%
  dplyr::left_join(df, by = "file_path") %>%
  left_join(conditionFriendlyNames()) %>%
  left_join(phenotypeFriendlyNames()) %>%
  dplyr::filter(parameter %in% c("eCLIP_both_ln","eCLIP_3end_ln"))


### Compare different quant methods
data = dplyr::filter(all_data, figure_name == "I", !is.na(phenotype)) %>%
  dplyr::filter(!(parameter %in% c("eCLIP_both_ln","eCLIP_3end_ln")))

methods_plot = ggplot(data, aes(x = estimate, y = phenotype, xmin = CI_lo, xmax = CI_hi)) + 
  geom_point() + facet_grid(parameter~.) +
  geom_errorbarh(aes(height = 0)) +
  theme_light() +
  xlab("Log2(enrichment)") +
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_cartesian(xlim = c(-.5,3.9))

ggsave("results/figures/fgwas_methods.pdf", plot = methods_plot, width = 5, height = 5)


#Compare different positioms
data = dplyr::filter(all_data, figure_name == "I", quant %in% 
                       c("txrevise_contained","txrevise_ends","txrevise_promoters","leafcutter")) %>%
  dplyr::filter(!(parameter %in% c("eCLIP_both_ln","eCLIP_3end_ln")))

methods_plot = ggplot(data, aes(x = estimate, y = quant, xmin = CI_lo, xmax = CI_hi)) + 
  geom_point() + facet_grid(parameter~.) +
  geom_errorbarh(aes(height = 0)) +
  theme_light() +
  xlab("Log2(enrichment)") +
  theme(strip.text.y = element_text(angle = 0)) + 
  coord_cartesian(xlim = c(-.5,3.9))
ggsave("results/figures/fgwas_position.pdf", plot = methods_plot, width = 5, height = 5)

