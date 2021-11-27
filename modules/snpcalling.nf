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
				[ -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log] && rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
				breseq --reference ${proteins} --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz
				"""
			else
				"""
				[ -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log] && rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
				breseq --reference ${proteins} --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz
				"""
		else
			if(params.mapping_breseq_p)
				"""
				[ -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log] && rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
				breseq --reference ${proteins} --num-processors $task.cpus --polymorphism-prediction --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz
				"""
			else
				"""
				[ -f ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log] && rm ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.log
				breseq --reference ${proteins} --num-processors $task.cpus --brief-html-output ${params.mapping_breseq_options} --output ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz
				"""
}

