import scipy.stats as st
import gzip

#st.norm.ppf(q, loc=0, scale=1)

qtltools_file = gzip.open('test.txt.gz', 'r')
fgwas_file = gzip.open('processed/annotations/fgwas/fgwas_annotations.txt.gz','r')
n_samples = '86'

#Make full header
header = "SNPID CHR POS Z F N SEGNUMBER"
fgwas_header = fgwas_file.readline().decode("utf8").rstrip().split("\t")
annotations = fgwas_header[5:]
new_header = " ".join([header] + annotations)
print(new_header)

#Initialize fgwas file
current_fgwas_pos = 0
current_fgwas_line = ""
current_fgwas_chr = "1"

for line in qtltools_file:
    line = line.decode("utf8").rstrip()
    fields = line.split("\t")
    gene_id = fields[0]
    snp_id = fields[7]
    snp_chr = fields[8]
    snp_pos = int(fields[9])
    pval = float(fields[11])
    effect = float(fields[12])
    af = '0'

    #Calculate z-score
    if(effect > 0):
        zscore = st.norm.ppf(pval/2, loc=0, scale=1)
    else:
        zscore = -st.norm.ppf(pval/2, loc=0, scale=1)

    #Iterate through the fgwas file
    if (current_fgwas_pos < snp_pos) & (snp_chr == current_fgwas_chr):
        while(current_fgwas_pos < snp_pos):
            current_fgwas_line = fgwas_file.readline().decode("utf8").rstrip()
            fields = current_fgwas_line.split("\t")
            current_fgwas_pos = int(fields[1])
            current_fgwas_chr = fields[0]
    elif snp_chr != current_fgwas_chr:
        while(snp_chr != current_fgwas_chr):
            current_fgwas_line = fgwas_file.readline().decode("utf8").rstrip()
            fields = current_fgwas_line.split("\t")
            current_fgwas_pos = int(fields[1])
            current_fgwas_chr = fields[0]
    else:
        while(current_fgwas_pos < snp_pos):
            current_fgwas_line = fgwas_file.readline().decode("utf8").rstrip()
            fields = current_fgwas_line.split("\t")
            current_fgwas_pos = int(fields[1])
            current_fgwas_chr = fields[0]

    #Extract relevant info from fgwas 
    if current_fgwas_pos == snp_pos:
        fgwas_fields = current_fgwas_line.split("\t")
        af = round(float(fgwas_fields[3]), 3)
        annots = fgwas_fields[5:]

    #Compile the final string
    string = " ".join([snp_id, snp_chr, str(snp_pos), str(zscore), str(af), n_samples, gene_id] + annots) 
    print(string)
