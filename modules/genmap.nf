process GENMAP {
	conda baseDir + '/env/snpless-genmap.yml'
	tag "GENMAP on ${reference}"
	cpus params.genmap_threads

	publishDir "${params.output}/GENMAP", mode: 'copy'

	input:
		path(reference)

	output:
		path(".")

	when:
		params.genmap || params.run_all

	script:
		if(params.run_genmap)
			"""
			genmap index -F ${reference} -I index/
			genmap map -K ${params.genmap_kmers} -E ${params.genmap_errors} -I index/ -O reference.genmap.K${params.genmap_kmers}.E${params.genmap_errors} ${params.genmap_outputs} -T $task.cpus
			"""
}
