#Construct new transcript annotations using reviseAnnotations
snakemake --cluster scripts/snakemake_submit.py -np -s reviseAnnotations.snakefile processed/annotations/reviseAnnotations/events/reviseAnnotations.transcript_metadata.txt --jobs 100 --configfile configs/reviseAnnotations_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s reviseAnnotations.snakefile processed/annotations/reviseAnnotations/events/reviseAnnotations.transcript_metadata.txt --jobs 100 --configfile configs/reviseAnnotations_config.GRCh38_87.yaml

#Construct promoter events
snakemake --cluster scripts/snakemake_submit.py -np -s txrevise_promoters.snakefile processed/annotations/txrevise_promoters/merged/txrevise_promoters.gff3 --jobs 100 --configfile configs/txrevise_promoters_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s txrevise_promoters.snakefile processed/annotations/txrevise/out.txt --jobs 100 --configfile configs/txrevise_promoters_config.yaml

#Quantify gene and transcript expression
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/salmonella/out.txt --jobs 100 --configfile configs/salmonella_trQTL_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s quantify_transcription.snakefile processed/acLDL/out.txt --jobs 100 --configfile configs/acLDL_trQTL_config.yaml

#Import transcript quants for acLDL data
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importTranscriptExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importReviseAnnotations.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importLeafCutter --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importLeafCutter.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importGeneExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/acLDL_importTxrevisePromoters.R"

#Import transcript quants for Salmonella data
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importTranscriptExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 28000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importReviseAnnotations.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importLeafCutter --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importLeafCutter.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 8000 --jobname importTranscriptExpression --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importGeneExpression.R"
echo "test" | python ~/software/utils/submitJobs.py --MEM 28000 --jobname importReviseAnnotations --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/salmonella_importTxrevisePromoters.R"

#Map trQTLs
snakemake --cluster scripts/snakemake_submit.py -np -s map_trQTLs.snakefile processed/salmonella/out.txt --jobs 20 --configfile configs/salmonella_trQTL_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s map_trQTLs.snakefile processed/acLDL/out.txt --jobs 1200 --configfile configs/acLDL_trQTL_config.yaml

#Run fgwas on all QTLs
snakemake --cluster scripts/snakemake_submit_UT.py -np -s run_fgwas.snakefile processed/salmonella/out.txt --jobs 30 --configfile configs/salmonella_trQTL_config.yaml
snakemake --cluster scripts/snakemake_submit_UT.py -np -s run_fgwas.snakefile processed/acLDL/out.txt --jobs 1200 --configfile configs/acLDL_trQTL_config.yaml

#Run coloc against all QTLs
snakemake --cluster scripts/snakemake_submit.py -np -s run_coloc.snakefile processed/salmonella/coloc_out.txt --jobs 500 --configfile configs/salmonella_trQTL_config.yaml
snakemake --cluster scripts/snakemake_submit.py -np -s run_coloc.snakefile processed/acLDL/coloc_out.txt --jobs 500 --configfile configs/acLDL_trQTL_config.yaml

#Convert revised GFFs into a single GRangesList object
echo "test" | python ~/software/utils/submitJobs.py --MEM 28000 --jobname importGFFs --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/munge/importRevisedGFFs.R"

#Estimate variance explained by main vs interaction effects for each condition"
echo "test" | python ~/software/utils/submitJobs.py --MEM 12000 --jobname varExp --ncores 1 --queue normal --command "/software/R-3.4.0/bin/Rscript analysis/trQTLs/estimate_condition_specificity.R"

#Run StringTie
snakemake --cluster scripts/snakemake_submit_UT.py -np -s run_StringTie.snakefile processed/salmonella/out.txt --jobs 100 --configfile configs/salmonella_trQTL_config.yaml


#Run fgwas to perform diease enrichments
snakemake --cluster scripts/snakemake_submit_UT.py -np -s run_fgwas_disease.snakefile processed/salmonella/out.txt --jobs 30 --configfile configs/salmonella_trQTL_config.yaml
