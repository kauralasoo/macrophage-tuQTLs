
#import variant information
GRCh37_information = importVariantInformation("results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz")
GRCh37_pos = dplyr::select(GRCh37_information, snp_id, chr, pos)

#Import lead qtls for each method
methods = c("featureCounts", "Ensembl_87", "leafcutter", "reviseAnnotations", 
            "txrevise_promoters", "txrevise_contained", "txrevise_ends")
method_list = setNames(as.list(methods), methods)
min_pvals = "processed/salmonella/fgwas/min_pvalues/"

qtltools_columns = c("phenotype_id","pheno_chr","pheno_start", "pheno_end",
                     "strand","n_snps", "distance", "snp_id", "snp_chr",
                     "snp_start", "snp_end", "p_nominal","beta", "is_lead")
qtl_list = purrr::map(method_list, ~readr::read_delim(paste0(min_pvals, .,"/naive.min_pvalues.txt.gz"), 
                                                      delim = "\t", col_names = qtltools_columns) %>% 
                        dplyr::filter(p_nominal < 1e-5))

#Makr QTL SNPs
qtl_snps = purrr::map(qtl_list, ~dplyr::select(., snp_id) %>% 
             dplyr::left_join(GRCh37_pos, by = "snp_id") %>% 
             dplyr::mutate(is_qtl = 1))
outdir = "processed/salmonella/garfield/annotations/"


makeAnnotationsByChr <- function(chrom, qtl_snps, outdir){
  #Run chr-by-chr
  chr_string = paste0("chr", chrom)
  print(chr_string)
  
  #Import Garfield coords
  variant_coords_dir = "~/projects/macrophage-gxe-study/databases/garfield-data/positions/"
  coords = readr::read_delim(file.path(variant_coords_dir, chr_string),
                             col_names = c("pos", "annot"), col_types = "ic", delim = " ") %>%
    dplyr::select(pos)
  
  #Make aannotations based on qtl_snps
  qtl_coords = purrr::map(qtl_snps, ~dplyr::filter(.,chr == chrom) %>% dplyr::select(pos, is_qtl))
  qtl_annots = purrr::map(qtl_coords, ~dplyr::left_join(coords, ., by = "pos") %>% 
                            dplyr::mutate(is_qtl = ifelse(is.na(is_qtl), 0, is_qtl)) %>% 
                            dplyr::select(is_qtl))
  all_annots = dplyr::bind_cols(coords, qtl_annots)
  all_annots$pos = paste0(all_annots$pos, " ")
  
  #Save annoations to disk
  write.table(all_annots, file.path(outdir, chr_string), sep = "", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
}


chromosomes = c(1:22) %>% as.character()
for (chr in chromosomes){
  makeAnnotationsByChr(chrom = chr, qtl_snps, outdir)
}

#Make a link file
link_file = data_frame(Index = c(0:6), annotation = methods, 
                       Celltype = "Macrophage", "Tissue" = "Blood", 
                       Type = rep("QTL",7), 
                       Category = "QTL")
write.table(link_file, paste0(outdir, "link_file.txt"), quote = FALSE,
            sep = " ", row.names = FALSE, col.names = TRUE)




