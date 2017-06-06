#Quantify gene and transcript expression
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/salmonella/out.txt --jobs 100 --configfile configs/salmonella_trQTL_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/acLDL/out.txt --jobs 100 --configfile configs/acLDL_trQTL_config.yaml

#Import transcript quants for acLDL data
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importTranscriptExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importReviseAnnotations.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importLeafCutter --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importLeafCutter.R"

#Import transcript quants for Salmonella data
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importTranscriptExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 28000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importReviseAnnotations.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importLeafCutter --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importLeafCutter.R"
