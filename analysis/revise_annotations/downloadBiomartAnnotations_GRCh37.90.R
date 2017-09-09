library("GenomicFeatures")
library("biomaRt")
library("dplyr")

#Make TranscriptDb object from biomart
txdb90 = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                             host="grch37.ensembl.org")
saveDb(txdb90, "../../annotations/GRCh37/genes/Ensembl_90/TranscriptDb_GRCh37_90.db")

#Downlaod transcript metadata from Ensembl
ensembl = useMart("ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)

#Define attributes to be downloaded from biomart
biomart_attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", 
                       "chromosome_name","strand", "transcript_start","transcript_end","ensembl_transcript_id", 
                       "transcript_version",
                       "transcript_gencode_basic","external_transcript_name", "source",
                       "transcript_length", "transcript_biotype", "ccds", "percentage_gene_gc_content")
refseq_attributes = c("ensembl_gene_id","ensembl_transcript_id","refseq_mrna")

#Download data for sample genes
genes = c("ENSG00000128604")
data = getBM(attributes = biomart_attributes, 
             filters = c("ensembl_gene_id"), 
             values = genes, 
             mart = ensembl)
data

data1 = getBM(attributes = refseq_attributes, 
              filters = c("ensembl_gene_id"), 
              values = genes, 
              mart = ensembl)
data1

#Download data for all transcripts
transcript_data = getBM(attributes = biomart_attributes, mart = ensembl)
refseq_data = getBM(attributes = refseq_attributes, mart = ensembl)
saveRDS(transcript_data, "../../annotations/GRCh37/genes/Ensembl_90/Homo_sapiens.GRCh37.90.transcript_data.rds")
saveRDS(refseq_data, "../../annotations/GRCh37/genes/Ensembl_90/Homo_sapiens.GRCh37.90.refseq_data.rds")

#Extract GENCODE basic transcript ids
gencode_basic = dplyr::tbl_df(transcript_data) %>% 
  dplyr::filter(transcript_gencode_basic == 1) %>% 
  dplyr::select(ensembl_gene_id, ensembl_transcript_id)
write.table(gencode_basic, "../../annotations/GRCh37/genes/Ensembl_90/Homo_sapiens.GRCh37.90.gencode_basic.txt", sep = "\t", row.names = FALSE, quote = FALSE)



