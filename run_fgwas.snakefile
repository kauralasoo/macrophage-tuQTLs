shell.prefix("set +o pipefail; ")

rule make_all:
	input:
		expand("processed/{{study}}/fgwas/output/{model}/{annot_type}/{condition}.params", annot_type = config["fgwas_phenotypes"], condition = config["conditions"], model = config["fgwas_model"]),
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Convert QTLtools output into format suitable for fgwas
rule make_fgwas_results:
	input:
		sorted = "processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz",
		perm = "processed/{study}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz"
	output:
		"processed/{study}/fgwas/input/{annot_type}/{condition}.fgwas_input.txt.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		source activate py3.6
		python scripts/qtltools_to_fgwas.py --qtltools {input.sorted} --annot {config[fgwas_annotations]} --N {config[sample_size]} --perm {input.perm} | gzip > {output}
		"""

rule sort_fgwas_results:
	input:
		"processed/{study}/fgwas/input/{annot_type}/{condition}.fgwas_input.txt.gz"
	output:
		"processed/{study}/fgwas/input/{annot_type}/{condition}.fgwas_input.sorted.txt.gz"
	resources:
		mem = 1000
	threads: 2
	shell:
		"""
		zcat {input} | head -n 1 | gzip > {output} && zcat {input} | tail -n +2 | sort -k7,7n -k2,2n -k3,3n | gzip >> {output}
		"""

rule run_fgwas:
	input:
		"processed/{study}/fgwas/input/{annot_type}/{condition}.fgwas_input.sorted.txt.gz"
	output:
		"processed/{study}/fgwas/output/{model}/{annot_type}/{condition}.params"
	params:
		out_prefix = "processed/{study}/fgwas/output/{model}/{annot_type}/{condition}"
	resources:
		mem = 56000
	threads: 1
	shell:
		"""
		module load gsl-1.16
		fgwas -i {input} -fine -o {params.out_prefix} -w {wildcards.model}
		"""
