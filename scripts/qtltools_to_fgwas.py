import scipy.stats as st
import gzip
import argparse

#st.norm.ppf(q, loc=0, scale=1)

parser = argparse.ArgumentParser(description = "Convert QTLtools output into format suitable for fgwas.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--qtltools", help = "Sorted QTLtools output file.")
parser.add_argument("--annot", help = "Variant annotations for fgwas (sorted by positition)")
parser.add_argument("--N", help = "QTL sample size.")
args = parser.parse_args()

#Set up input files
qtltools_file = gzip.open(args.qtltools, 'r')
fgwas_file = gzip.open(args.annot,'r')
n_samples = args.N

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

#Gene dictionary
gene_dict = dict()
gene_count = 0

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

    #Convert gene id to SEGNUMBER
    if gene_id in gene_dict:
        gene_number = gene_dict[gene_id]
    else:
        gene_count = gene_count + 1  
        gene_number = gene_count
        gene_dict[gene_id] = gene_number

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
    string = " ".join([snp_id, snp_chr, str(snp_pos), str(zscore), str(af), n_samples, str(gene_number)] + annots) 
    print(string)
