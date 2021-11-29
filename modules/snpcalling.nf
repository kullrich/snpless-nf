process FREEBAYESBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BRESEQ.freebayes.vcf"

	when:
		(params.snpcalling && params.run_freebayes) || params.run_all

	script:
		"""
		freebayes -f ${reference} ${params.snpcalling_freebayes_options} -b ${mappingPath}/*.bam > BRESEQ.freebayes.vcf
		"""
}

process FREEBAYESMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.freebayes.vcf"

	when:
		(params.snpcalling && params.run_freebayes) || params.run_all

	script:
		"""
		freebayes -f ${reference} ${params.snpcalling_freebayes_options} -b ${mappingPath}/*.bam > MINIMAP2.freebayes.vcf
		"""
}

process FREEBAYESBWA {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.freebayes.vcf"

	when:
		(params.snpcalling && params.run_freebayes) || params.run_all

	script:
		"""
		freebayes -f ${reference} ${params.snpcalling_freebayes_options} -b ${mappingPath}/*.bam > BWA.freebayes.vcf
		"""
}
