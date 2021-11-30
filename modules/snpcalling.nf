process FREEBAYESBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/FREEBAYES/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.freebayes.vcf"

	when:
		(params.snpcalling && params.run_freebayes) || params.run_all

	script:
		"""
		freebayes -f ${mappingPath}/reference.fasta ${params.snpcalling_freebayes_options} -b ${mappingPath}/*.bam > BRESEQ.freebayes.vcf
		"""
}

process FREEBAYESMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/FREEBAYES/MINIMAP2", mode: 'symlink'

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

	publishDir "${params.output}/SNPCALLING/FREEBAYES/BWA", mode: 'symlink'

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

process BCFTOOLSBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-bcftools.yml'
	tag "BCFTOOLS on ${mappingPath}"
	cpus params.bcftools_threads

	publishDir "${params.output}/SNPCALLING/BCFTOOLS/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.bcftools.vcf"

	when:
		(params.snpcalling && params.run_bcftools) || params.run_all

	script:
		"""
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${mappingPath}/reference.fasta ${mappingPath}/*.bam | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > BRESEQ.bcftools.vcf
		"""
}

process BCFTOOLSMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-bcftools.yml'
	tag "BCFTOOLS on ${mappingPath}"
	cpus params.bcftools_threads

	publishDir "${params.output}/SNPCALLING/BCFTOOLS/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.bcftools.vcf"

	when:
		(params.snpcalling && params.run_bcftools) || params.run_all

	script:
		"""
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${reference} ${mappingPath}/*.bam | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > MINIMAP2.bcftools.vcf
		"""
}

process BCFTOOLSBWA {
	conda baseDir + '/env/snpless-snpcalling-bcftools.yml'
	tag "BCFTOOLS on ${mappingPath}"
	cpus params.bcftools_threads

	publishDir "${params.output}/SNPCALLING/BCFTOOLS/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.bcftools.vcf"

	when:
		(params.snpcalling && params.run_bcftools) || params.run_all

	script:
		"""
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${reference} ${mappingPath}/*.bam | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > BWA.bcftools.vcf
		"""
}

process GDCOMPARE {
	conda baseDir + '/env/snpless-mapping-breseq.yml'
	tag "GDCOMPARE on ${mappingPath}"
	cpus params.gdtools_threads

	publishDir "${params.output}/SNPCALLING/GDCOMPARE/", mode: 'symlink'

	input:
		val(gds)
		path(mappingPath)

	output:
		path "BRESEQ.annotate.*"

	when:
		(params.snpcalling && params.run_gdcompare) || params.run_all

	script:
		"""
		gdtools ANNOTATE -o BRESEQ.annotate.html -r ${mappingPath}/reference.fasta ${mappingPath}/*.gd
		gdtools ANNOTATE -f TSV -o BRESEQ.annotate.tsv -r ${mappingPath}/reference.fasta ${mappingPath}/*.gd
		"""
}
