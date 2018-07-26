#Convert gff3 files generated by reviseAnnotations into fasta sequence
rule convert_gff3_to_fasta:
	input:
		"processed/annotations/gff/{annotation}.gff3"
	output:
		"processed/annotations/fasta/{annotation}.fa"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		module load cufflinks2.2
		module load python-2.7.13
		gffread -w {output} -g {config[reference_genome]} {input}
		"""


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
rule quant_salmon:
	input:
		fq1 = "processed/{study}/shuffled/{sample}_1.fq.gz",
		fq2 = "processed/{study}/shuffled/{sample}_2.fq.gz",
		salmon_index = "processed/annotations/salmon_index/{annotation}"
	output:
		"processed/{study}/salmon/{annotation}/{sample}/quant.sf"
	params:
		out_prefix = "processed/{study}/salmon/{annotation}/{sample}"
	resources:
		mem = 2000
	threads: 2	
	shell:
		"salmon --no-version-check quant --libType {config[libType]} "
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

rule make_fastq:
	input:
		fa1 = "processed/{study}/fastq/{sample}_1.fasta.gz",
		fa2 = "processed/{study}/fastq/{sample}_2.fasta.gz"
	output:
		fq1 = "processed/{study}/fq/{sample}_1.fq.gz",
		fq2 = "processed/{study}/fq/{sample}_2.fq.gz"
	resources:
		mem = 100
	threads: 2
	shell:
		"""
		module load perl-5.22.0
		zcat {input.fa1} | perl scripts/fasta_to_fastq.pl - | gzip > {output.fq1}
		zcat {input.fa2} | perl scripts/fasta_to_fastq.pl - | gzip > {output.fq2}
		"""

rule shuffle_fastq:
	input:
		fq1 = "processed/{study}/fq/{sample}_1.fq.gz",
		fq2 = "processed/{study}/fq/{sample}_2.fq.gz"
	output:
		fq1 = "processed/{study}/shuffled/{sample}_1.fq.gz",
		fq2 = "processed/{study}/shuffled/{sample}_2.fq.gz"
	resources:
		mem = 6000
	threads: 1
	shell:
		"""
		~/anaconda3/envs/py3.6/bin/seqkit shuffle {input.fq1} -s 100 | gzip > {output.fq1}
		~/anaconda3/envs/py3.6/bin/seqkit shuffle {input.fq2} -s 100 | gzip > {output.fq2}
		"""

rule quantify_whippet:
	input:
		fq1 = "processed/{study}/shuffled/{sample}_1.fq.gz",
		fq2 = "processed/{study}/shuffled/{sample}_2.fq.gz"
	output:
		"processed/{study}/whippet/{sample}.psi.gz"
	params:
		out = "processed/{study}/whippet/{sample}"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		module load julia-0.6.0
		julia ~/.julia/v0.6/Whippet/bin/whippet-quant.jl {input.fq1} {input.fq2} -o {params.out} -x {config[whippet_index]}
		"""

rule hisat2_align:
	input:
		fq1 = "processed/{study}/shuffled/{sample}_1.fq.gz",
		fq2 = "processed/{study}/shuffled/{sample}_2.fq.gz"
	output:
		bam = "processed/{study}/hisat2/{sample}.bam"
	resources:
		mem = 8000
	threads: 2
	shell:
		"""
		module load samtools-1.6
		hisat2 -p {threads} -x {config[hisat2_index]} {config[hisat2_flags]} -1 {input.fq1} -2 {input.fq2} | samtools view -Sb > {output.bam}
		"""

rule quantify_dexseq:
	input:
		bam = "processed/{study}/hisat2/{sample}.bam"
	output:
		counts = "processed/{study}/DEXseq/{sample}.counts"
	resources:
		mem = 2000
	threads: 1
	shell:
		"""
		source activate py2.7
		python scripts/dexseq_count.py {config[dexseq_gtf]} {input.bam} {output.counts} -p yes -s no -f bam -r name
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{study}}/matrices/{annotation}.salmon_txrevise.rds", annotation=config["annotations"]),
		expand("processed/{{study}}/whippet/{sample}.psi.gz", sample=config["samples"]),
		expand("processed/{{study}}/hisat2/{sample}.bam", sample=config["samples"]),
		expand("processed/{{study}}/DEXseq/{sample}.counts", sample=config["samples"])
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"