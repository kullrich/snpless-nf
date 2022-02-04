process FASTQC {
	conda baseDir + '/env/snpless-qc-fastqc.yml'
	tag "FASTQC on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.fastqc_threads

	publishDir "${params.output}/QC/FASTQC", mode: 'symlink'

	input:
		tuple val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reads1f), path(reads2f), path(trim_adapater_file)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.qc && params.run_fastqc) || params.run_all

	script:
		if (reads2 == "-")
			if (params.qc_fastqc_quiet)
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
			else
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
		else
			if (params.qc_fastqc_quiet)
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
			else
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.fastqc.log
					"""
}

process TRIM {
	conda baseDir + '/env/snpless-qc-trim.yml'
	tag "TRIM on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.trim_threads

	publishDir "${params.output}/QC/TRIM", mode: 'symlink'

	input:
		tuple val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(reads1f), path(reads2f), path(trim_adapter_file)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.qc && params.run_trim) || params.run_all

	script:
		if (reads2 == "-")
			if (params.qc_fastqc_quiet)
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic SE -threads $task.cpus -quiet $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic SE -threads $task.cpus -quiet $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					"""
			else
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic SE -threads $task.cpus $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic SE -threads $task.cpus $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					"""
		else
			if (params.qc_fastqc_quiet)
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic PE -threads $task.cpus -quiet $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic PE -threads $task.cpus -quiet $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq
					"""
			else
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic PE -threads $task.cpus $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}
					trimmomatic PE -threads $task.cpus $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.trim.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE2.fq
					"""
}

process PEAR {
	conda baseDir + '/env/snpless-qc-pear.yml'
	tag "PEAR on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.pear_threads

	publishDir "${params.output}/QC/PEAR", mode: 'symlink'

	input:
		tuple path(trimPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)


	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.qc && params.run_pear) || params.run_all

	script:
		if (reads2 == "-")
			"""
			"""
		else
			"""
			pear --threads $task.cpus ${params.qc_pear_options} -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE1.fq.gz -r ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.PE2.fq.gz -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.pear.log
			mv ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.unassembled.forward.fastq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq
			mv ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.unassembled.reverse.fastq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.assembled.fastq
			gunzip -c ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.assembled.fastq ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz > ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.assembled.fastq.gz
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.discarded.fastq
			"""
}
