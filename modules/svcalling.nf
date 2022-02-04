process PINDELBRESEQ {
	conda baseDir + '/env/snpless-svcalling-pindel.yml'
	tag "PINDELBRESEQ on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.pindel_threads

	publishDir "${params.output}/SVCALLING/PINDEL/BRESEQ", mode: 'symlink'

	input:
		tuple path(breseqPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: pindel_breseq
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.pindel*", emit: pindel_files

	when:
		(params.svcalling && params.run_pindel && params.run_breseq && !params.skip_breseq && !params.skip_pindel) || params.run_all

	script:
		"""
		samtools view ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.sorted.bam | sam2pindel - ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.pindel.txt ${params.svcalling_pindel_sam2pindel_options}
		samtools faidx ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/reference.fasta
		pindel -T $task.cpus -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/reference.fasta -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.pindel.txt ${params.svcalling_pindel_options} -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.pindel
		"""
}

process PINDELMINIMAP2 {
	conda baseDir + '/env/snpless-svcalling-pindel.yml'
	tag "PINDELMINIMAP2 on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.pindel_threads

	publishDir "${params.output}/SVCALLING/PINDEL/MINIMAP2", mode: 'symlink'

	input:
		tuple path(minimap2Path), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)
		path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: pindel_minimap2
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.pindel*", emit: pindel_files

	when:
		(params.svcalling && params.run_pindel && params.run_minimap2 && !params.skip_minimap2 && !params.skip_pindel) || params.run_all

	script:
		"""
		samtools view ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam | sam2pindel - ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.pindel.txt ${params.svcalling_pindel_sam2pindel_options}
		samtools faidx ${reference}
		pindel -T $task.cpus -f ${reference} -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.pindel.txt ${params.svcalling_pindel_options} -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.pindel
		"""
}

process PINDELBWA {
	conda baseDir + '/env/snpless-svcalling-pindel.yml'
	tag "PINDELBWA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.pindel_threads

	publishDir "${params.output}/SVCALLING/PINDEL/BWA", mode: 'symlink'

	input:
		tuple path(bwaPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)
		path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: pindel_bwa
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.pindel*", emit: pindel_files

	when:
		(params.svcalling && params.run_pindel && params.run_bwa && !params.skip_bwa && !params.skip_pindel) || params.run_all

	script:
		"""
		samtools view ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam | sam2pindel - ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.pindel.txt ${params.svcalling_pindel_sam2pindel_options}
		samtools faidx ${reference}
		pindel -T $task.cpus -f ${reference} -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.pindel.txt ${params.svcalling_pindel_options} -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.pindel
		"""
}

process GRIDSSBRESEQ {
	conda baseDir + '/env/snpless-svcalling-gridss.yml'
	tag "VARSCANBRESEQ on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.varscan.*.vcf"

	when:
		(params.snpcalling && params.run_varscan) || params.run_all

	script:
		"""
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${mappingPath}/reference.fasta ${mappingPath}/*.bam | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > BRESEQ.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${mappingPath}/reference.fasta ${mappingPath}/*.bam | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > BRESEQ.varscan.indel.vcf
		"""
}

process GRIDSSMINIMAP2 {
	conda baseDir + '/env/snpless-svcalling-gridss.yml'
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
		(params.snpcalling && params.run_varscan) || params.run_all

	script:
		"""
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > MINIMAP2.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > MINIMAP2.varscan.indel.vcf
		"""
}

process GRIDSSBWA {
	conda baseDir + '/env/snpless-svcalling-gridss.yml'
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
		(params.snpcalling && params.run_varscan) || params.run_all

	script:
		"""
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > BWA.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} ${mappingPath}/*.bam | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > BWA.varscan.indel.vcf
		"""
}
