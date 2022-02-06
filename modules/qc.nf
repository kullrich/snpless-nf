process FASTQC {
	conda baseDir + '/env/snpless-qc-fastqc.yml'
	tag "FASTQC on ${sampleLong}"
	cpus params.fastqc_threads

	publishDir "${params.output}/QC/FASTQC", mode: 'symlink'

	input:
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reads1f), path(reads2f), path(trim_adapater_file)

	output:
		path("${sampleLong}"), emit: fastqcDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: fastqcSamples


	when:
		(params.qc && params.run_fastqc) || params.run_all

	script:
		if (reads2 == "-")
			if (params.qc_fastqc_quiet)
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q --nogroup -o ${sampleLong} $reads1f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q -o ${sampleLong} $reads1f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
			else
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers --nogroup -o ${sampleLong} $reads1f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -o ${sampleLong} $reads1f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
		else
			if (params.qc_fastqc_quiet)
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q --nogroup -o ${sampleLong} $reads1f $reads2f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q -o ${sampleLong} $reads1f $reads2f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
			else
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers --nogroup -o ${sampleLong} $reads1f $reads2f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleLong}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -o ${sampleLong} $reads1f $reads2f &> ${sampleLong}/${sampleLong}.fastqc.log
					"""
}

process TRIM {
	conda baseDir + '/env/snpless-qc-trim.yml'
	tag "TRIM on ${sampleLong}"
	cpus params.trim_threads

	publishDir "${params.output}/QC/TRIM", mode: 'symlink'

	input:
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reads1f), path(reads2f), path(trim_adapter_file)

	output:
		path("${sampleLong}"), emit: trimDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: trimSamples

	when:
		(params.qc && params.run_trim) || params.run_all

	script:
		if (reads2 == "-")
			if (params.qc_fastqc_quiet)
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleLong}
					trimmomatic SE -threads $task.cpus -quiet $reads1f ${sampleLong}/${sampleLong}.SE.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					"""
				else
					"""
					mkdir -p ${sampleLong}
					trimmomatic SE -threads $task.cpus -quiet $reads1f ${sampleLong}/${sampleLong}.SE.fq ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					"""
			else
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleLong}
					trimmomatic SE -threads $task.cpus $reads1f ${sampleLong}/${sampleLong}.SE.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					"""
				else
					"""
					mkdir -p ${sampleLong}
					trimmomatic SE -threads $task.cpus $reads1f ${sampleLong}/${sampleLong}.SE.fq ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					"""
		else
			if (params.qc_fastqc_quiet)
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleLong}
					trimmomatic PE -threads $task.cpus -quiet $reads1f $reads2f ${sampleLong}/${sampleLong}.PE1.fq ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.PE2.fq ${sampleLong}/${sampleLong}.SE2.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.PE1.fq
					bgzip ${sampleLong}/${sampleLong}.PE2.fq 
					cat ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq > ${sampleLong}/${sampleLong}.SE.fq
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					rm ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq
					"""
				else
					"""
					mkdir -p ${sampleLong}
					trimmomatic PE -threads $task.cpus -quiet $reads1f $reads2f ${sampleLong}/${sampleLong}.PE1.fq ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.PE2.fq ${sampleLong}/${sampleLong}.SE2.fq ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.PE1.fq
					bgzip ${sampleLong}/${sampleLong}.PE2.fq 
					cat ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq > ${sampleLong}/${sampleLong}.SE.fq
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					rm ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq
					"""
			else
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleLong}
					trimmomatic PE -threads $task.cpus $reads1f $reads2f ${sampleLong}/${sampleLong}.PE1.fq ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.PE2.fq ${sampleLong}/${sampleLong}.SE2.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.PE1.fq
					bgzip ${sampleLong}/${sampleLong}.PE2.fq 
					cat ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq > ${sampleLong}/${sampleLong}.SE.fq
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					rm ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq
					"""
				else
					"""
					mkdir -p ${sampleLong}
					trimmomatic PE -threads $task.cpus $reads1f $reads2f ${sampleLong}/${sampleLong}.PE1.fq ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.PE2.fq ${sampleLong}/${sampleLong}.SE2.fq ${params.qc_trim_options} &> ${sampleLong}/${sampleLong}.trim.log
					bgzip ${sampleLong}/${sampleLong}.PE1.fq
					bgzip ${sampleLong}/${sampleLong}.PE2.fq 
					cat ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq > ${sampleLong}/${sampleLong}.SE.fq
					bgzip ${sampleLong}/${sampleLong}.SE.fq
					rm ${sampleLong}/${sampleLong}.SE1.fq ${sampleLong}/${sampleLong}.SE2.fq
					"""
}

process PEAR {
	conda baseDir + '/env/snpless-qc-pear.yml'
	tag "PEAR on ${sampleLong}"
	cpus params.pear_threads

	publishDir "${params.output}/QC/PEAR", mode: 'symlink'

	input:
		path(trimDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		path("${sampleLong}"), emit: pearDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: pearSamples

	when:
		(params.qc && params.run_pear) || params.run_all

	script:
		if (reads2 == "-")
			"""
			"""
		else
			"""
			pear --threads $task.cpus ${params.qc_pear_options} -f ${sampleLong}/${sampleLong}.PE1.fq.gz -r ${sampleLong}/${sampleLong}.PE2.fq.gz -o ${sampleLong}/${sampleLong} &> ${sampleLong}/${sampleLong}.pear.log
			mv ${sampleLong}/${sampleLong}.unassembled.forward.fastq ${sampleLong}/${sampleLong}.MEPE1.fq
			bgzip ${sampleLong}/${sampleLong}.MEPE1.fq
			mv ${sampleLong}/${sampleLong}.unassembled.reverse.fastq ${sampleLong}/${sampleLong}.MEPE2.fq
			bgzip ${sampleLong}/${sampleLong}.MEPE2.fq
			bgzip ${sampleLong}/${sampleLong}.assembled.fastq
			gunzip -c ${sampleLong}/${sampleLong}.assembled.fastq ${sampleLong}/${sampleLong}.SE.fq.gz > ${sampleLong}/${sampleLong}.MESE.fq
			bgzip ${sampleLong}/${sampleLong}.MESE.fq
			rm ${sampleLong}/${sampleLong}.assembled.fastq.gz
			rm ${sampleLong}/${sampleLong}.discarded.fastq
			"""
}
