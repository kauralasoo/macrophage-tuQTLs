import gzip
import argparse
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description = "Convert uniform GWAS summary stats into format suitable for fgwas.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--gwas", help = "GWAS summary stats file.")
parser.add_argument("--annot", help = "Variant annotations for fgwas (sorted by positition)")
args = parser.parse_args()

#Open files
gwas_file = gzip.open(args.gwas, 'r')
fgwas_file = gzip.open(args.annot,'r')

#Make new header 
header = gwas_file.readline().decode("utf8").rstrip()
fgwas_header = fgwas_file.readline().decode("utf8").rstrip().split("\t")
header = "SNPID CHR POS F Z N SE"
annotations = fgwas_header[5:]
new_header = " ".join([header] + annotations)
print(new_header)

#Initialize fgwas file
current_fgwas_pos = 0
current_fgwas_line = ""
current_fgwas_chr = "1"

for line in gwas_file:
    line = line.decode("utf8").rstrip()
    fields = line.split()
    snp_id = fields[0]
    snp_chr = fields[1]
    snp_pos = fields[2]
    freq = 1
    zscore = fields[10]
    n_samples = 10000
    se = fields[9]
    if se == "NA":
        se = '1'

    #print(fields)
    string = " ".join([snp_id, snp_chr, snp_pos, str(freq), zscore, str(n_samples), se])
    print(string)