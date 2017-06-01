#Quantify gene and transcript expression
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/salmonella/out.txt --jobs 100 --configfile configs/salmonella_trQTL_config.yaml
