

#Import transcript metadata
tx_meta = readRDS("results/simulations/transcript_meta.rds")
gene_transcript_map = dplyr::transmute(tx_meta, gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id) %>%
  dplyr::filter(transcript_id %in% all_tx)

#Import truncated annotations
exons = readRDS("results/simulations/both_transcripts.rds")

#Import QTLs
vcf_file = readRDS("results/genotypes/salmonella/imputed.86_samples.sorted.filtered.named.rds")
salmonella_qtls = readRDS("results/trQTLs/salmonella_trQTL_min_pvalues.rds")

#Import QTL effects
extended_by_truncation = saveRDS("results/simulations/sim_one_effects.rds")

#Import lead variants
lead_variants = dplyr::filter(salmonella_qtls$Ensembl_87$naive, group_id %in% tx_meta$ensembl_gene_id) %>%
  dplyr::transmute(gene_id = group_id, snp_id)

track_data = dplyr::data_frame(sample_id = colnames(vcf_file$genotypes), genotype_id = colnames(vcf_file$genotypes), 
                               bigWig = list.files("processed/sim_one/bigwig/", full.names = T), scaling_factor = 1, track_id = "simulation")



#First example
selected_gene = "ENSG00000187147"
selected_tx = dplyr::filter(gene_transcript_map, gene_id == "ENSG00000187147")$transcript_id
selected_snp_id = "rs12125215"

#Make coverage plot
geno_track_data = wiggleplotrGenotypeColourGroup(track_data, selected_snp_id, vcf_file$genotypes, -1)


coverage_plot = plotCoverage(exons = exons[selected_tx], cdss = NULL, 
                             track_data = geno_track_data, rescale_introns = T,
                             fill_palette = getGenotypePalette(), 
                             plot_fraction = 0.3, heights = c(0.5,0.5), 
                             return_subplots_list = F, coverage_type = "line", transcript_label = TRUE)
ggsave("results/figures/simulated_coverage.pdf", plot = coverage_plot, width = 7, height = 4)
coverage_plot


#Make boxplots as well
event_ids = dplyr::filter(extended_by_truncation, gene_id == "ENSG00000187147")$transcript_id

# Import tx usage
#repeat for extended transcripts
upstream_extended = readRDS("processed/sim_one/matrices/txrevise_upstream_both.salmon_txrevise.rds")
downstream_extended = readRDS("processed/sim_one/matrices/txrevise_downstream_both.salmon_txrevise.rds")

#Construct event annotations
upstream_events = data_frame(transcript_id = rownames(upstream_extended$abundance)) %>%
  tidyr::separate(transcript_id, c("gene_id", "suffix"), sep = "\\.grp_1", remove = FALSE) %>%
  dplyr::select(-suffix)

downstream_events = data_frame(transcript_id = rownames(downstream_extended$abundance)) %>%
  tidyr::separate(transcript_id, c("gene_id", "suffix"), sep = "\\.grp_1", remove = FALSE) %>%
  dplyr::select(-suffix)

#Calculate event ratios
tu_up_original = calculateTranscriptRatios(upstream_extended$abundance, upstream_events)
tu_down_original = calculateTranscriptRatios(downstream_extended$abundance, downstream_events)
event_usage = rbind(tu_up_original, tu_down_original)

#Construct joint metadata
txrevise_metadata = rbind(upstream_events, downstream_events) %>% dplyr::arrange(gene_id) %>%
  dplyr::left_join(dplyr::transmute(tx_meta, gene_id = ensembl_gene_id, gene_name = external_gene_name) %>% dplyr::distinct()) %>% 
  dplyr::transmute(gene_id = transcript_id, gene_name)

#Upstream event boxplot
up_event_id = "ENSG00000187147.grp_1.upstream.ENST00000484745"
plot_data = constructQtlPlotDataFrame(up_event_id, selected_snp_id, 
                                      event_usage, vcf_file$genotypes, sample_metadata, txrevise_metadata)
up_usage = ggplot(plot_data, aes(x = genotype_value, y = norm_exp)) + geom_boxplot() + geom_point() +
  theme_light() +
  ylab("ENST00000484745 promoter usage") + 
  xlab("genotype") +
  coord_cartesian(ylim = c(0,1))
ggsave("results/figures/simulated_up_usage.pdf", plot = up_usage, width = 3, height = 3)


#3'end event boxplot
down_event_id = "ENSG00000187147.grp_1.downstream.ENST00000484745"
plot_data = constructQtlPlotDataFrame(down_event_id, selected_snp_id, 
                                      event_usage, vcf_file$genotypes, sample_metadata, txrevise_metadata)
down_usage = ggplot(plot_data, aes(x = genotype_value, y = norm_exp)) + geom_boxplot() + geom_point() +
  theme_light() +
  ylab("ENST00000484745 3'end usage") + 
  xlab("genotype") +
  coord_cartesian(ylim = c(0,1))
ggsave("results/figures/simulated_down_usage.pdf", plot = down_usage, width = 3, height = 3)



