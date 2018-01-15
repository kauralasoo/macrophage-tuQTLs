import scipy.stats as st
import gzip

#st.norm.ppf(q, loc=0, scale=1)

qtltools_file = gzip.open('test.txt.gz', 'r')
fgwas_file = gzip.open('processed/annotations/fgwas/fgwas_annotations.txt.gz','r')

for line in qtltools_file:
    line = line.decode("utf8").rstrip()
    fields = line.split("\t")
    gene_id = fields[0]
    snp_id = fields[7]
    snp_chr = fields[8]
    snp_pos = fields[9]
    string = " ".join([snp_id, snp_chr, snp_pos, "0", "0", "0", gene_id]) 
    print(string)
