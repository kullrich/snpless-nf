process BRESEQ {
	conda baseDir + '/env/snpless-mapping-breseq.yml'
	tag "BRESEQ on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.breseq_threads

	publishDir "${params.output}/MAPPING/BRESEQ", mode: 'copy'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), path(proteins)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_breseq) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_breseq_p)
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz
				"""
			else
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz
				"""
		else
			if(params.mapping_breseq_p)
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz
				"""
			else
				"""
				breseq --reference ${proteins} --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz
				"""
}

process MINIMAP2 {
	conda baseDir + '/env/snpless-mapping-minimap2.yml'
	tag "MINMAP2 on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.minimap2_threads

	publishDir "${params.output}/MAPPING/MINIMAP2", mode: 'copy'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_minimap2) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_minimap2_samblaster)
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.bam
				"""
			else
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.bam
				"""
		else
			if(params.mapping_minimap2_samblaster)
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster -r | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.paired.bam
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.single.bam
				"""
			else
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.paired.bam
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.single.bam
				"""
}

process BWA {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "BWA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.bwa_threads

	publishDir "${params.output}/MAPPING/BWA", mode: 'copy'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), path(reference)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.assembly && params.run_prokka) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_bwa_samblaster)
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r -M | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.bam
				"""
			else
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.bam
				"""
		else
			if(params.mapping_bwa_samblaster)
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster -r -M | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.paired.bam
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r -M | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.single.bam
				"""
			else
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.paired.bam
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz -R '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.single.bam
				"""
}

process POSTBRESEQ {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "POSTBRESEQ on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.breseq_threads

	publishDir "${params.output}/MAPPING/BRESEQ", mode: 'copy'

	input:
		tuple path(breseqPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	output:
		tuple path("BRESEQOUT"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_breseq) || params.run_all

	script:
		"""
		mkdir -p BRESEQOUT
		cp ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/output.vcf BRESEQOUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.vcf
		cp ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/output.gd BRESEQOUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.gd
		samtools addreplacerg -r '@RG\\tID:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tSM:${sampleId}_${sampleReplicate}_${sampleTimepoint}\\tLB:lib1\\tPL:illumina\\tPU:unit1' -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}/data/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bam BRESEQOUT/reference.bam
		samtools index BRESEQOUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bam
		samtools coverage ${params.mapping_breseq_coverage} BRESEQOUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bam > BRESEQOUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.mean.coverage
		"""
}

process POSTMINIMAP2 {
	conda baseDir + '/env/snpless-mapping-minimap2.yml'
	tag "POSTMINIMAP2 on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.minimap2_threads

	publishDir "${params.output}/MAPPING/MINIMAP2", mode: 'copy'

	input:
		tuple path(minimap2Path), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	output:
		tuple path("MINIMAP2OUT"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_minimap2) || params.run_all

	script:
		if (reads2 == "-")
			"""
			mkdir -p MINIMAP2OUT
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.bam > MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam
			samtools index MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam
			samtools coverage ${params.mapping_minimap2_coverage} MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam > MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.mean.coverage
			"""
		else
			"""
			mkdir -p MINIMAP2OUT
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.paired.bam > MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.paired.sorted.bam
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.single.bam > MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.single.sorted.bam
			samtools merge -@ $task.cpus MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.paired.sorted.bam MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.single.sorted.bam
			rm MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.paired.sorted.bam 
			rm MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.single.sorted.bam
			samtools index MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam
			samtools coverage ${params.mapping_minimap2_coverage} MINIMAP2OUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.minimap2.mean.coverage
			"""
}

process POSTBWA {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "POSTBWA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.bwa_threads

	publishDir "${params.output}/MAPPING/BWAOUT/${sampleId}_${sampleReplicate}_${sampleTimepoint}", mode: 'copy'

	input:
		tuple path(bwaPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	output:
		tuple path("."), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.mapping && params.run_bwa) || params.run_all

	script:
		if (reads2 == "-")
			"""
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam
			samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam
			samtools coverage ${params.mapping_bwa_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.mean.coverage
			"""
		else
			"""
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.paired.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.paired.sorted.bam
			samtools sort -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.single.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.single.sorted.bam
			samtools merge -@ $task.cpus ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.paired.sorted.bam ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.single.sorted.bam
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.paired.sorted.bam 
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.single.sorted.bam
			samtools index ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam
			samtools coverage ${params.mapping_bwa_coverage} ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.sorted.bam > ${sampleId}_${sampleReplicate}_${sampleTimepoint}.bwa.mean.coverage
			"""
}
