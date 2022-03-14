process BRESEQ {
	conda baseDir + '/env/snpless-mapping-breseq.yml'
	tag "BRESEQ on ${sampleLong}"
	cpus params.breseq_threads

	publishDir "${params.output}/MAPPING/BRESEQ", mode: 'symlink'

	input:
		path(pearDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reference), path(gff3)

	output:
		path("${sampleLong}"), emit: breseqDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: breseqSamples

	when:
		(params.mapping && params.run_breseq && !params.skip_breseq) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_breseq_p)
				"""
				(cat ${gff3}; echo "##FASTA"; cat ${reference}) > reference.gff3
				breseq --reference reference.gff3 --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleLong} ${sampleLong}/${sampleLong}.SE.fq.gz
				"""
			else
				"""
				(cat ${gff3}; echo "##FASTA"; cat ${reference}) > reference.gff3
				breseq --reference reference.gff3 --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleLong} ${sampleLong}/${sampleLong}.SE.fq.gz
				"""
		else
			if(params.mapping_breseq_p)
				"""
				(cat ${gff3}; echo "##FASTA"; cat ${reference}) > reference.gff3
				breseq --reference reference.gff3 --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleLong} ${sampleLong}/${sampleLong}.MEPE1.fq.gz ${sampleLong}/${sampleLong}.MEPE2.fq.gz ${sampleLong}/${sampleLong}.MESE.fq.gz
				"""
			else
				"""
				(cat ${gff3}; echo "##FASTA"; cat ${reference}) > reference.gff3
				breseq --reference reference.gff3 --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleLong} ${sampleLong}/${sampleLong}.MEPE1.fq.gz ${sampleLong}/${sampleLong}.MEPE2.fq.gz ${sampleLong}/${sampleLong}.MESE.fq.gz
				"""
}

process MINIMAP2 {
	conda baseDir + '/env/snpless-mapping-minimap2.yml'
	tag "MINMAP2 on ${sampleLong}"
	cpus params.minimap2_threads

	publishDir "${params.output}/MAPPING/MINIMAP2", mode: 'symlink'

	input:
		path(pearDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reference)

	output:
		path("${sampleLong}"), emit: minimap2Dir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: minimap2Samples

	when:
		(params.mapping && params.run_minimap2 && !params.skip_minimap2) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_minimap2_samblaster)
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleLong}/${sampleLong}.SE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r | samtools view -b - > ${sampleLong}/${sampleLong}.minimap2.bam
				"""
			else
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleLong}/${sampleLong}.SE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleLong}/${sampleLong}.minimap2.bam
				"""
		else
			if(params.mapping_minimap2_samblaster)
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleLong}/${sampleLong}.MEPE1.fq.gz ${sampleLong}/${sampleLong}.MEPE2.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster -r | samtools view -b - > ${sampleLong}/${sampleLong}.minimap2.paired.bam
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleLong}/${sampleLong}.MESE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r | samtools view -b - > ${sampleLong}/${sampleLong}.minimap2.single.bam
				"""
			else
				"""
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleLong}/${sampleLong}.MEPE1.fq.gz ${sampleLong}/${sampleLong}.MEPE2.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleLong}/${sampleLong}.minimap2.paired.bam
				minimap2 -t $task.cpus ${params.mapping_minimap2_options} ${reference} ${sampleLong}/${sampleLong}.MESE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleLong}/${sampleLong}.minimap2.single.bam
				"""
}

process BWA {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "BWA on ${sampleLong}"
	cpus params.bwa_threads

	publishDir "${params.output}/MAPPING/BWA", mode: 'symlink'

	input:
		path(pearDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reference)

	output:
		path("${sampleLong}"), emit: bwaDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: bwaSamples

	when:
		(params.mapping && params.run_bwa && !params.skip_bwa) || params.run_all

	script:
		if (reads2 == "-")
			if(params.mapping_bwa_samblaster)
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleLong}/${sampleLong}.SE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r -M | samtools view -b - > ${sampleLong}/${sampleLong}.bwa.bam
				"""
			else
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleLong}/${sampleLong}.SE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleLong}/${sampleLong}.bwa.bam
				"""
		else
			if(params.mapping_bwa_samblaster)
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleLong}/${sampleLong}.MEPE1.fq.gz ${sampleLong}/${sampleLong}.MEPE2.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster -r -M | samtools view -b - > ${sampleLong}/${sampleLong}.bwa.paired.bam
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleLong}/${sampleLong}.MESE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samblaster --ignoreUnmated -r -M | samtools view -b - > ${sampleLong}/${sampleLong}.bwa.single.bam
				"""
			else
				"""
				bwa index ${reference}
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleLong}/${sampleLong}.MEPE1.fq.gz ${sampleLong}/${sampleLong}.MEPE2.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleLong}/${sampleLong}.bwa.paired.bam
				bwa mem -t $task.cpus ${params.mapping_bwa_options} ${reference} ${sampleLong}/${sampleLong}.MESE.fq.gz -R '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' | samtools view -b - > ${sampleLong}/${sampleLong}.bwa.single.bam
				"""
}

process POSTBRESEQ {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "POSTBRESEQ on ${sampleLong}"
	cpus params.breseq_threads

	publishDir "${params.output}/MAPPING/BRESEQ", mode: 'copy'

	input:
		path(breseqDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		path("${sampleLong}"), emit: postbreseqDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: postbreseqSamples
		path "${sampleLong}/data/${sampleLong}.breseq.mean.coverage", emit: breseq_mean_coverage
		path "${sampleLong}/data/${sampleLong}.breseq.sorted.bam", emit: breseq_bam
		path "${sampleLong}/data/${sampleLong}.breseq.sorted.bam.bai", emit: breseq_bam_index
		path "${sampleLong}/data/reference.fasta", emit: breseq_bam_reference
		path "${sampleLong}/data/reference.gff3", emit: breseq_bam_reference_gff3
		path "${sampleLong}/data/${sampleLong}.vcf", emit: breseq_vcf
		path "${sampleLong}/data/${sampleLong}.GT.vcf", emit: breseq_gtvcf
		path "${sampleLong}/data/${sampleLong}.gd", emit: breseq_gd

	when:
		(params.mapping && params.run_breseq && !params.skip_breseq) || params.run_all

	script:
		"""
		cp ${sampleLong}/data/output.vcf ${sampleLong}/data/${sampleLong}.vcf
		cp ${sampleLong}/data/output.gd ${sampleLong}/data/${sampleLong}.gd
		samtools addreplacerg -r '@RG\\tID:${sampleLong}\\tSM:${sampleLong}\\tLB:lib1\\tPL:illumina\\tPU:unit1' -o ${sampleLong}/data/${sampleLong}.breseq.bam ${sampleLong}/data/reference.bam
		samtools sort -@ $task.cpus ${sampleLong}/data/${sampleLong}.breseq.bam  > ${sampleLong}/data/${sampleLong}.breseq.sorted.bam 
		samtools index ${sampleLong}/data/${sampleLong}.breseq.sorted.bam
		samtools coverage ${params.mapping_breseq_coverage} ${sampleLong}/data/${sampleLong}.breseq.sorted.bam > ${sampleLong}/data/${sampleLong}.breseq.mean.coverage
		breseqVCFaddGT.py ${sampleLong}/data/${sampleLong}.vcf ${sampleLong} > ${sampleLong}/data/${sampleLong}.GT.vcf
		"""
}

process POSTMINIMAP2 {
	conda baseDir + '/env/snpless-mapping-minimap2.yml'
	tag "POSTMINIMAP2 on ${sampleLong}"
	cpus params.minimap2_threads

	publishDir "${params.output}/MAPPING/MINIMAP2", mode: 'copy'

	input:
		path(minimap2Dir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		path("${sampleLong}"), emit: postminimap2Dir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: postminimap2Samples
		path "${sampleLong}/${sampleLong}.minimap2.mean.coverage", emit: minimap2_mean_coverage
		path "${sampleLong}/${sampleLong}.minimap2.sorted.bam", emit: minimap2_bam
		path "${sampleLong}/${sampleLong}.minimap2.sorted.bam.bai", emit: minimap2_bam_index

	when:
		(params.mapping && params.run_minimap2 && !params.skip_minimap2) || params.run_all

	script:
		if (reads2 == "-")
			"""
			samtools sort -@ $task.cpus ${sampleLong}/${sampleLong}.minimap2.bam > ${sampleLong}/${sampleLong}.minimap2.sorted.bam
			samtools index ${sampleLong}/${sampleLong}.minimap2.sorted.bam
			samtools coverage ${params.mapping_minimap2_coverage} ${sampleLong}/${sampleLong}.minimap2.sorted.bam > ${sampleLong}/${sampleLong}.minimap2.mean.coverage
			"""
		else
			"""
			samtools sort -@ $task.cpus ${sampleLong}/${sampleLong}.minimap2.paired.bam > ${sampleLong}/${sampleLong}.minimap2.paired.sorted.bam
			samtools sort -@ $task.cpus ${sampleLong}/${sampleLong}.minimap2.single.bam > ${sampleLong}/${sampleLong}.minimap2.single.sorted.bam
			samtools merge -@ $task.cpus ${sampleLong}/${sampleLong}.minimap2.sorted.bam ${sampleLong}/${sampleLong}.minimap2.paired.sorted.bam ${sampleLong}/${sampleLong}.minimap2.single.sorted.bam
			rm ${sampleLong}/${sampleLong}.minimap2.paired.sorted.bam 
			rm ${sampleLong}/${sampleLong}.minimap2.single.sorted.bam
			samtools index ${sampleLong}/${sampleLong}.minimap2.sorted.bam
			samtools coverage ${params.mapping_minimap2_coverage} ${sampleLong}/${sampleLong}.minimap2.sorted.bam > ${sampleLong}/${sampleLong}.minimap2.mean.coverage
			"""
}

process POSTBWA {
	conda baseDir + '/env/snpless-mapping-bwa.yml'
	tag "POSTBWA on ${sampleLong}"
	cpus params.bwa_threads

	publishDir "${params.output}/MAPPING/BWA", mode: 'copy'

	input:
		path(bwaDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		path("${sampleLong}"), emit: postbwaDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: postbwaSamples
		path "${sampleLong}/${sampleLong}.bwa.mean.coverage", emit: bwa_mean_coverage
		path "${sampleLong}/${sampleLong}.bwa.sorted.bam", emit: bwa_bam
		path "${sampleLong}/${sampleLong}.bwa.sorted.bam.bai", emit: bwa_bam_index

	when:
		(params.mapping && params.run_bwa && !params.skip_bwa) || params.run_all

	script:
		if (reads2 == "-")
			"""
			samtools sort -@ $task.cpus ${sampleLong}/${sampleLong}.bwa.bam > ${sampleLong}/${sampleLong}.bwa.sorted.bam
			samtools index ${sampleLong}/${sampleLong}.bwa.sorted.bam
			samtools coverage ${params.mapping_bwa_coverage} ${sampleLong}/${sampleLong}.bwa.sorted.bam > ${sampleLong}/${sampleLong}.bwa.mean.coverage
			"""
		else
			"""
			samtools sort -@ $task.cpus ${sampleLong}/${sampleLong}.bwa.paired.bam > ${sampleLong}/${sampleLong}.bwa.paired.sorted.bam
			samtools sort -@ $task.cpus ${sampleLong}/${sampleLong}.bwa.single.bam > ${sampleLong}/${sampleLong}.bwa.single.sorted.bam
			samtools merge -@ $task.cpus ${sampleLong}/${sampleLong}.bwa.sorted.bam ${sampleLong}/${sampleLong}.bwa.paired.sorted.bam ${sampleLong}/${sampleLong}.bwa.single.sorted.bam
			rm ${sampleLong}/${sampleLong}.bwa.paired.sorted.bam 
			rm ${sampleLong}/${sampleLong}.bwa.single.sorted.bam
			samtools index ${sampleLong}/${sampleLong}.bwa.sorted.bam
			samtools coverage ${params.mapping_bwa_coverage} ${sampleLong}/${sampleLong}.bwa.sorted.bam > ${sampleLong}/${sampleLong}.bwa.mean.coverage
			"""
}
