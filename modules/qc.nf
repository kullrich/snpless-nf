process FASTQC {
	conda 'env/snpless-qc-fastqc.yml'
	tag "FASTQC on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.fastqc_threads

	publishDir "${params.output}/QC/FASTQC", mode: 'copy'

	input:
		tuple val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), path(reads1f), path(reads2f), path(trim_adapater_file)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.qc && params.run_fastqc) || params.run_all

	script:
		if (reads2 == "-")
			if (params.qc_fastqc_quiet)
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
			else
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
		else
			if (params.qc_fastqc_quiet)
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -q -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
			else
				if (params.qc_fastqc_nogroup)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers --nogroup -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					fastqc -t $task.cpus -k $params.qc_fastqc_kmers -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} $reads1f $reads2f &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					"""
}

process TRIM {
	conda 'env/snpless-qc-trim.yml'
	tag "TRIM on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.trim_threads

	publishDir "${params.output}/QC/TRIM", mode: 'copy'

	input:
		tuple val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), path(reads1f), path(reads2f), path(trim_adapter_file)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.qc && params.run_trim) || params.run_all

	script:
		if (reads2 == "-")
			if (params.qc_fastqc_quiet)
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic SE -threads $task.cpus -quiet $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic SE -threads $task.cpus -quiet $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					"""
			else
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic SE -threads $task.cpus $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic SE -threads $task.cpus $reads1f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					"""
		else
			if (params.qc_fastqc_quiet)
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic PE -threads $task.cpus -quiet $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic PE -threads $task.cpus -quiet $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq
					"""
			else
				if(params.qc_trim_use_adapter)
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic PE -threads $task.cpus $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq ILLUMINACLIP:${trim_adapter_file}:${params.qc_clip_options} ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq
					"""
				else
					"""
					mkdir -p ${sampleId}_${sampleReplicate}_${sampleTimepoint}
					trimmomatic PE -threads $task.cpus $reads1f $reads2f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq ${params.qc_trim_options} &> ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq 
					cat ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq
					rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE1.fq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE2.fq
					"""
}

process PEAR {
	conda 'env/snpless-qc-pear.yml'
	tag "PEAR on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.pear_threads

	publishDir "${params.output}/QC/PEAR", mode: 'copy'

	input:
		tuple path(trimPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)


	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.qc && params.run_pear) || params.run_all

	script:
		if (reads2 == "-")
			"""
			"""
		else
			"""
			pear --threads $task.cpus ${params.qc_pear_options} -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq.gz -r ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq.gz -o ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint} > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
			mv ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.unassembled.forward.fastq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq
			mv ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.unassembled.reverse.fastq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.assembled.fastq
			gunzip -c ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.assembled.fastq ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz > ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq
			bgzip ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE1.fq.gz
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.PE2.fq.gz
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.assembled.fastq.gz
			rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.discarded.fastq
			"""
}
