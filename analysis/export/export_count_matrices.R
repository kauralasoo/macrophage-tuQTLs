library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("devtools")
library("SummarizedExperiment")

#AcLDL read counts
acLDL_fc = readRDS("results/SummarizedExperiments/acLDL_featureCounts.rds")
gz1 = gzfile("results/tables/featureCounts/acLDL_count_matrix.txt.gz", "w")
write.table(assays(acLDL_fc)$counts, gz1, sep = "\t", quote = FALSE, row.names = T, col.names = T)
close(gz1)

gz1 = gzfile("results/tables/featureCounts/acLDL_cqn_matrix.txt.gz", "w")
write.table(assays(acLDL_fc)$cqn, gz1, sep = "\t", quote = FALSE, row.names = T, col.names = T)
close(gz1)

gz1 = gzfile("results/tables/featureCounts/acLDL_gene_metadata.txt.gz", "w")
write.table(rowData(acLDL_fc), gz1, sep = "\t", quote = FALSE, row.names = F, col.names = T)
close(gz1)

gz1 = gzfile("results/tables/featureCounts/acLDL_sample_metadata.txt.gz", "w")
write.table(colData(acLDL_fc), gz1, sep = "\t", quote = FALSE, row.names = F, col.names = T)
close(gz1)

#Salmonella read counts
acLDL_fc = readRDS("results/SummarizedExperiments/salmonella_featureCounts.rds")
gz1 = gzfile("results/tables/featureCounts/salmonella_count_matrix.txt.gz", "w")
write.table(assays(acLDL_fc)$counts, gz1, sep = "\t", quote = FALSE, row.names = T, col.names = T)
close(gz1)

gz1 = gzfile("results/tables/featureCounts/salmonella_cqn_matrix.txt.gz", "w")
write.table(assays(acLDL_fc)$cqn, gz1, sep = "\t", quote = FALSE, row.names = T, col.names = T)
close(gz1)

gz1 = gzfile("results/tables/featureCounts/salmonella_gene_metadata.txt.gz", "w")
write.table(rowData(acLDL_fc), gz1, sep = "\t", quote = FALSE, row.names = F, col.names = T)
close(gz1)

gz1 = gzfile("results/tables/featureCounts/salmonella_sample_metadata.txt.gz", "w")
write.table(colData(acLDL_fc), gz1, sep = "\t", quote = FALSE, row.names = F, col.names = T)
close(gz1)