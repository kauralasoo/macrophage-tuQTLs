#Quantify gene and transcript expression
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/salmonella/out.txt --jobs 100 --configfile configs/salmonella_trQTL_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/salmonella/out.txt --jobs 100 --configfile configs/acLDL_trQTL_config.yaml

#Import transcript quants
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importTranscriptExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript macrophage-gxe-study/acLDL/munge/acLDL_importReviseAnnotations.R"
