rule run_all:
	input:
		expand("processed/{{study}}/coloc/{gwas}.{phenotype}.{coloc_window}.txt", gwas = config["gwas_traits"], phenotype = config["coloc_phenotypes"], coloc_window = config["coloc_window"])
	output:
		"processed/{study}/coloc_out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Run coloc accross inflammatory traits
rule run_coloc:
	input:
		"processed/{study}/qtltools/output/{phenotype}/sorted"
	output:
		"processed/{study}/coloc/{gwas}.{phenotype}.{coloc_window}.txt"
	params:
		outdir = "processed/{study}/coloc",
		qtl_dir = "processed/{study}/qtltools/output/",
		phenotype = "{phenotype}",
		gwas = "{gwas}",
		coloc_window = "{coloc_window}"

	resources:
		mem = 12000
	threads: 1
	shell:
		"/software/R-3.4.0/bin/Rscript analysis/coloc/GWAS_run_coloc.R "
		"--phenotype {wildcards.phenotype} --window {wildcards.coloc_window} "
		"--gwas {wildcards.gwas} --dir {config[gwas_dir]} --outdir {params.outdir} "
		"--qtl {params.qtl_dir} --samplesizes {config[sample_sizes]} "
		"--gwasvarinfo {config[gwas_var_info]} --qtlvarinfo {config[qtl_var_info]} "
		"--gwaslist {config[gwas_list]}"

