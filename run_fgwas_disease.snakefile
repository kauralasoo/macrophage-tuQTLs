rule make_all:
	input:
		expand("processed/{{study}}/fgwas/min_pvalues/{annot_type}/{condition}.min_pvalues.txt.gz", annot_type = config["fgwas_phenotypes"], condition = config["conditions"]),
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"


#Extract minimal p-values for each position from QTLtools
rule extract_min_pvals:
	input:
		sorted = "processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz",
		perm = "processed/{study}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz"
	output:
		"processed/{study}/fgwas/min_pvalues/{annot_type}/{condition}.min_pvalues.txt.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		source activate py3.6
		python scripts/qtltools_keep_min_pvalue.py --qtltools {input.sorted} --perm {input.perm} | gzip > {output}
		"""
