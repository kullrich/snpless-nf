process BRESEQ {
	conda baseDir + '/env/snpless-mapping-breseq.yml'
	tag "BRESEQ on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.breseq_threads

	publishDir "${params.output}/MAPPING/BRESEQ", mode: 'symlink'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(proteins)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_breseq && !params.skip_breseq) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_breseq_p)
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz
				"""
			else
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz
				"""
		else
			if(params.mapping_breseq_p)
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz
				"""
			else
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz
				"""
}

process MINIMAP2 {
	conda baseDir + '/env/snpless-mapping-minimap2.yml'
	tag "MINMAP2 on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.minimap2_threads

	publishDir "${params.output}/MAPPING/MINIMAP2", mode: 'symlink'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_minimap2 && !params.skip_minimap2) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_minimap2_samblaster)
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.bam
				"""
			else
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.bam
				"""
		else
			if(params.mapping_minimap2_samblaster)
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster -r | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.paired.bam
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.single.bam
				"""
			else
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.paired.bam
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.single.bam
				"""
}

process BWA {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "BWA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.bwa_threads

	publishDir "${params.output}/MAPPING/BWA", mode: 'symlink'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_bwa && !params.skip_bwa) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_bwa_samblaster)
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r -M | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.bam
				"""
			else
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.bam
				"""
		else
			if(params.mapping_bwa_samblaster)
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster -r -M | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.paired.bam
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r -M | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.single.bam
				"""
			else
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.paired.bam
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.single.bam
				"""
}

process POSTBRESEQ {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "POSTBRESEQ on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.breseq_threads

	publishDir "${params.output}/MAPPING/BRESEQ", mode: 'copy'

	input:
		tuple path(breseqPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: postbreseq
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.mean.coverage", emit: breseq_mean_coverage
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.sorted.bam", emit: breseq_bam
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.sorted.bam.bai", emit: breseq_bam_index
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/reference.fasta", emit: breseq_bam_reference
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/reference.gff3", emit: breseq_bam_reference_gff3
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.vcf", emit: breseq_vcf
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.gd", emit: breseq_gd

	when:
		(params.mapping && params.run_breseq && !params.skip_breseq) || params.run_all

	script:
		"""
		cp ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/output.vcf ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.vcf
		cp ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/output.gd ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.gd
		samtools addreplacerg -r '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}\\tLB:lib1\\tPL:illumina\\tPU:unit1' -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/reference.bam
		samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.bam  > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.sorted.bam 
		samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.sorted.bam
		samtools coverage ${params.mapping_breseq_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.breseq.mean.coverage
		"""
}

process POSTMINIMAP2 {
	conda baseDir + '/env/snpless-mapping-minimap2.yml'
	tag "POSTMINIMAP2 on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.minimap2_threads

	publishDir "${params.output}/MAPPING/MINIMAP2", mode: 'copy'

	input:
		tuple path(minimap2Path), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: postminimap2
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.mean.coverage", emit: minimap2_mean_coverage
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam", emit: minimap2_bam
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam.bai", emit: minimap2_bam_index

	when:
		(params.mapping && params.run_minimap2 && !params.skip_minimap2) || params.run_all

	script:
		if (reads2 == "-")
			"""
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam
			samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam
			samtools coverage ${params.mapping_minimap2_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.mean.coverage
			"""
		else
			"""
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.paired.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.paired.sorted.bam
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.single.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.single.sorted.bam
			samtools merge -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.paired.sorted.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.single.sorted.bam
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.paired.sorted.bam 
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.single.sorted.bam
			samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam
			samtools coverage ${params.mapping_minimap2_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.minimap2.mean.coverage
			"""
}

process POSTBWA {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "POSTBWA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.bwa_threads

	publishDir "${params.output}/MAPPING/BWA", mode: 'copy'

	input:
		tuple path(bwaPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: postbwa
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.mean.coverage", emit: bwa_mean_coverage
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam", emit: bwa_bam
		path "${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam.bai", emit: bwa_bam_index

	when:
		(params.mapping && params.run_bwa && !params.skip_bwa) || params.run_all

	script:
		if (reads2 == "-")
			"""
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam
			samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam
			samtools coverage ${params.mapping_bwa_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.mean.coverage
			"""
		else
			"""
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.paired.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.paired.sorted.bam
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.single.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.single.sorted.bam
			samtools merge -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.paired.sorted.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.single.sorted.bam
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.paired.sorted.bam 
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.single.sorted.bam
			samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam
			samtools coverage ${params.mapping_bwa_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.bwa.mean.coverage
			"""
}
