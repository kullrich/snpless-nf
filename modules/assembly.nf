process UNICYCLER {
	conda 'env/snpless-assembly-unicycler.yml'
	tag "UNICYCLER on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.unicycler_threads

	publishDir "${params.output}/ASSEMBLY/UNICYCLER", mode: 'copy'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.assembly && params.run_unicycler) || params.run_all

	script:
		if (reads2 == "-")
			"""
			unicycler -s ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.SE.fq.gz -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} -t $task.cpus --mode ${params.assembly_unicycler_mode}
			mv ${sampleId}_${sampleReplicate}_${sampleTimepoint}/assembly.fasta ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.assembly.fasta
			"""
		else
			"""
			unicycler -1 ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE1.fq.gz -2 ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MEPE2.fq.gz -s ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.MESE.fq.gz -o ${sampleId}_${sampleReplicate}_${sampleTimepoint} -t $task.cpus --mode ${params.assembly_unicycler_mode}
			mv ${sampleId}_${sampleReplicate}_${sampleTimepoint}/assembly.fasta ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.assembly.fasta
			"""
}

process PROKKA {
	conda 'env/snpless-assembly-prokka.yml'
	tag "PROKKA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}"
	cpus params.prokka_threads

	publishDir "${params.output}/ASSEMBLY/PROKKA", mode: 'copy'

	input:
		tuple path(unicyclerPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2), path(proteins)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(reads1), val(reads2)

	when:
		(params.assembly && params.run_prokka) || params.run_all

	script:
		"""
		prokka --proteins ${proteins} --outdir ${sampleId}_${sampleReplicate}_${sampleTimepoint}/prokka --prefix ${sampleId}_${sampleReplicate}_${sampleTimepoint} ${sampleId}_${sampleReplicate}_${sampleTimepoint}/${sampleId}_${sampleReplicate}_${sampleTimepoint}.assembly.fasta
		"""
}
