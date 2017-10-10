#Are upstream events more likely to condition-specific
inter_df = dplyr::bind_rows(salmonella_df)
inter_genes = dplyr::filter(inter_df, interaction_fraction > 0.5, p_fdr < 0.1) %>%
  dplyr::filter(!is.na(naive_tpm)) %>%
  dplyr::mutate(not_expressed = ifelse(naive_tpm > 2, FALSE, TRUE)) %>%
  dplyr::mutate(shared_eQTL = ifelse(R2 > 0.8, TRUE, FALSE)) %>%
  dplyr::mutate(is_extreme_DE = ifelse(log2FoldChange < 5 | is.na(log2FoldChange), FALSE, TRUE)) %>%
  dplyr::mutate(expression_only = not_expressed | is_extreme_DE) %>%
  dplyr::mutate(all_technical = shared_eQTL | not_expressed | is_extreme_DE) %>%
  dplyr::mutate()

revised_qtls = dplyr::filter(salmonella_df, quant == "reviseAnnotations") %>%
  tidyr::separate(phenotype_id, c("ensembl_gene_id","group_id", "position", "transcript_id"), sep = "\\.", remove = FALSE)
revised_trqtls = dplyr::filter(inter_genes, quant == "reviseAnnotations") %>%
  tidyr::separate(phenotype_id, c("ensembl_gene_id","group_id", "position", "transcript_id"), sep = "\\.", remove = FALSE)

table(revised_qtls$position)
table(revised_trqtls$position)