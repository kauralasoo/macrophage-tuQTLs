#Build salmon indexes for fasta files
rule construct_salmon_index:
	input:
		"processed/annotations/fasta/{annotation}.fa"
	output:
		"processed/annotations/salmon_index/{annotation}"
	resources:
		mem = 10000
	threads: 1
	shell:
		"salmon -no-version-check index -t {input} -i {output}"

#Quantify gene expression using full Ensembl annotations
rule reviseAnnotation_quant_salmon:
	input:
		fq1 = "processed/{study}/fastq/{sample}_1.fastq.gz",
		fq2 = "processed/{study}/fastq/{sample}_2.fastq.gz",
		salmon_index = "processed/annotations/salmon_index/{annotation}"
	output:
		"processed/{study}/salmon/{annotation}/{sample}/quant.sf"
	params:
		out_prefix = "processed/{study}/salmon/{annotation}/{sample}"
	resources:
		mem = 10000
	threads: 8	
	shell:
		"salmon --no-version-check quant --seqBias --gcBias --libType {config[libType]} "
		"--index {input.salmon_index} -1 {input.fq1} -2 {input.fq2} -p {threads} "
		"-o {params.out_prefix}"

#Merge Salmon results
rule merge_salmon:
	input:
		expand("processed/{{study}}/salmon/{{annotation}}/{sample}/quant.sf", sample=config["samples"])
	output:
		"processed/{study}/matrices/{annotation}.salmon_txrevise.rds"
	params:
		sample_ids = ','.join(config["samples"]),
		dir = "processed/{study}/salmon/{annotation}"
	threads: 1
	resources:
		mem = 12000
	shell:
		"""
		module load R-3.4.1
		Rscript scripts/merge_Salmon.R -s {params.sample_ids} -d {params.dir} -o {output}
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{study}}/matrices/{annotation}.salmon_txrevise.rds", annotation=config["annotations"]),
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"