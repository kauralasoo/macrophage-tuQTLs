import scipy.stats as st
import gzip
import argparse
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description = "Keep one minimal p-value per position to make fgwas annotations.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--qtltools", help = "Sorted QTLtools output file.")
parser.add_argument("--perm", help = "QTLtools output from the permutation run.")
args = parser.parse_args()

#Set up input files
qtltools_file = gzip.open(args.qtltools, 'r')
perm_file = gzip.open(args.perm,'r')

#Make a directory of phenotypes to be included in the fgwas output
phenotype_dict = dict()
for line in perm_file:
    line = line.decode("utf8").rstrip()
    fields = line.split()
    phenotype_id = fields[5]
    phenotype_dict[phenotype_id] = 1
perm_file.close()

last_line = ""
last_snp = "" 
last_p = float(1)

#Iterate over the qtltools file
for line in qtltools_file:
    line = line.decode("utf8").rstrip()
    fields = line.split("\t")
    gene_id = fields[0]
    snp_id = fields[7]
    snp_chr = fields[8]
    snp_pos = int(fields[9])
    pval = float(fields[11])
    effect = float(fields[12])

    #Exit loop if gene_id not in phentype dict
    if not(gene_id in phenotype_dict):
        continue

    #Initialise
    if last_snp == "":
        last_snp = snp_id
        last_p = pval
        last_line = line

    #Keep minimal p-value per position
    if last_snp == snp_id:
        if pval < last_p:
            last_line = line
    else:
        print(last_line)
        last_snp = snp_id
        last_p = pval
        last_line = line

#Print the last line that was not printed yet
print(list_line)




