library("polyester")
library("Biostrings")

fold_changes = matrix(c(4,4,rep(1,904),1,1,0.25,0.25,rep(1,902)), nrow=906)
head(fold_changes)

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet("results/simulations/original_transcripts.fa")

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:10]
writeXStringSet(small_fasta, 'chr22_small.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)


# simulation call:
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
                    num_reps=c(3,3), fold_changes=fold_changes[877:886,], outdir='simulated_reads') 


