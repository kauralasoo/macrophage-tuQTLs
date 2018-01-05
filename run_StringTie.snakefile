#Run StringTie on each bam file
rule run_StringTie:
	input:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		"processed/{study}/StringTie/{sample}.gtf"
	params:
		local_tmp = "/tmp/" + uuid.uuid4().hex + "/"
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		module load samtools-1.6
		mkdir {params.local_tmp}
		cp {input} {params.local_tmp}/{wildcards.sample}.bam
        samtools view -h {params.local_tmp}/{wildcards.sample}.bam | gawk -v strType=2 -f tagXSstrandedData.awk | samtools view -bS - > {params.local_tmp}/{wildcards.sample}.XS.bam
		stringtie {params.local_tmp}/{wildcards.sample}.XS.bam --rf -o {output} -p {threads}
		rm -r {params.local_tmp}
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{study}/StringTie/{sample}.gtf", study = config["study"], sample=config["samples"]),
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"
