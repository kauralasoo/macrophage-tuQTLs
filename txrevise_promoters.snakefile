#Run reviseannotations for a single batch
rule run_txrevise:
	input:
		txdb = config["txdb"]
	output:
		gff = "processed/annotations/txrevise/{position}/txrevise_{position}.batch_{batch}_{n_batches}.gff3",
	params:
		out_prefix = "processed/annotations/txrevise/{position}/",
		chunk = "'{batch} {n_batches}'",
		pos = "{position}"
	threads: 1
	resources:
		mem = 2000
	shell:
		"/software/R-3.4.0/bin/Rscript analysis/revise_annotations/constructTruePromoterEvents.R -t {input.txdb} "
		"-b {params[chunk]} -o {params[out_prefix]} -p {params.pos}"

#Iterate over batches
rule merge_txrevise_batches:
	input:
		gff = expand("processed/annotations/txrevise/{{position}}/txrevise_{{position}}.batch_{batch}_{n_batches}.gff3",
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"]),
	output:
		gff = "processed/annotations/txrevise/{position}/merged/txrevise_{position}.gff3",
	resources:
		mem = 100
	threads: 1
	shell:
		'cat {input.gff} | grep -v "^#" > {output.gff}'

rule make_all:
	input: 
		expand("processed/annotations/txrevise/{position}/merged/txrevise_{position}.gff3", position = config["position"])
	output:
		"processed/annotations/txrevise/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"



