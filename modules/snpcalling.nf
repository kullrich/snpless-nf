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
		(params.snpcalling && params.run_freebayes && params.run_breseq && !params.skip_breseq && !params.skip_freebayes) || params.run_all

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
		(params.snpcalling && params.run_freebayes && params.run_minimap2 && !params.skip_minimap2 && !params.skip_freebayes) || params.run_all

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
		(params.snpcalling && params.run_freebayes && params.run_bwa && !params.skip_bwa && !params.skip_freebayes) || params.run_all

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
		(params.snpcalling && params.run_bcftools && params.run_breseq && !params.skip_breseq && !params.skip_bcftools) || params.run_all

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
		(params.snpcalling && params.run_bcftools && params.run_minimap2 && !params.skip_minimap2 && !params.skip_bcftools) || params.run_all

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
		(params.snpcalling && params.run_bcftools && params.run_bwa && !params.skip_bwa && !params.skip_bcftools) || params.run_all

	script:
		"""
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${reference} ${mappingPath}/*.bam | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > BWA.bcftools.vcf
		"""
}

process LOFREQBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-lofreq.yml'
	tag "LOFREQBRESEQ on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.lofreq_threads

	publishDir "${params.output}/SNPCALLING/LOFREQ/BRESEQ", mode: 'symlink'

	input:
		tuple path(breseqPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), emit: lofreq_breseq
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.breseq.lofreq.vcf", emit: lofreq_vcf

	when:
		(params.snpcalling && params.run_lofreq && params.run_breseq && !params.skip_breseq && !params.skip_lofreq) || params.run_all

	script:
		"""
		lofreq faidx ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/reference.fasta
		lofreq indelqual --dindel -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/reference.fasta ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}.breseq.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}.breseq.sorted.indel.bam
		lofreq index ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}.breseq.sorted.indel.bam
		lofreq call-parallel --pp-threads $task.cpus ${params.snpcalling_lofreq_options} -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/reference.fasta -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.breseq.lofreq.vcf ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}.breseq.sorted.indel.bam
		"""
}

process LOFREQMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-lofreq.yml'
	tag "LOFREQMINIMAP2 on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.lofreq_threads

	publishDir "${params.output}/SNPCALLING/LOFREQ/MINIMAP2", mode: 'symlink'

	input:
		tuple path(minimap2Path), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)
		path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), emit: lofreq_minimap2
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.lofreq.vcf", emit: lofreq_vcf

	when:
		(params.snpcalling && params.run_lofreq && params.run_minimap2 && !params.skip_minimap2 && !params.skip_lofreq) || params.run_all

	script:
		"""
		lofreq faidx ${reference}
		lofreq indelqual --dindel -f ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.indel.bam
		lofreq index ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.indel.bam
		lofreq call-parallel --pp-threads $task.cpus ${params.snpcalling_lofreq_options} -f ${reference} -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.lofreq.vcf ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.indel.bam
		"""
}

process LOFREQBWA {
	conda baseDir + '/env/snpless-snpcalling-lofreq.yml'
	tag "LOFREQBWA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.lofreq_threads

	publishDir "${params.output}/SNPCALLING/LOFREQ/BWA", mode: 'symlink'

	input:
		tuple path(bwaPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)
		path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), emit: lofreq_bwa
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.lofreq.vcf", emit: lofreq_vcf

	when:
		(params.snpcalling && params.run_lofreq && params.run_bwa && !params.skip_bwa && !params.skip_lofreq) || params.run_all

	script:
		"""
		lofreq faidx ${reference}
		lofreq indelqual --dindel -f ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.indel.bam
		lofreq index ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.indel.bam
		lofreq call-parallel --pp-threads $task.cpus ${params.snpcalling_lofreq_options} -f ${reference} -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.lofreq.vcf ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.indel.bam
		"""
}

process VARSCANBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "VARSCANBRESEQ on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.varscan.*.vcf"

	when:
		(params.snpcalling && params.run_varscan && params.run_breseq && !params.skip_breseq && !params.skip_varscan) || params.run_all

	script:
		"""
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${mappingPath}/reference.fasta ${mappingPath}/*.bam | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > BRESEQ.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${mappingPath}/reference.fasta ${mappingPath}/*.bam | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > BRESEQ.varscan.indel.vcf
		"""
}

process VARSCANMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "VARSCANMINIMAP2 on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.varscan.*.vcf"

	when:
		(params.snpcalling && params.run_varscan && params.run_minimap2 && !params.skip_minimap2 && !params.skip_varscan) || params.run_all

	script:
		"""
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > MINIMAP2.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > MINIMAP2.varscan.indel.vcf
		"""
}

process VARSCANBWA {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "VARSCANBWA on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.varscan.*.vcf"

	when:
		(params.snpcalling && params.run_varscan && params.run_bwa && !params.skip_bwa && !params.skip_varscan) || params.run_all

	script:
		"""
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > BWA.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > BWA.varscan.indel.vcf
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
		path(proteins)

	output:
		path "BRESEQ.annotate.*"

	when:
		(params.snpcalling && params.run_gdcompare && params.run_breseq && !params.skip_breseq && !params.skip_gdcompare) || params.run_all

	script:
		"""
		gdtools ANNOTATE -o BRESEQ.annotate.html -r ${proteins} ${mappingPath}/*.gd
		gdtools ANNOTATE -f TSV -o BRESEQ.annotate.tsv -r ${proteins} ${mappingPath}/*.gd
		"""
}
