process UNICYCLER {
	conda baseDir + '/env/snpless-assembly-unicycler.yml'
	tag "UNICYCLER on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.unicycler_threads

	publishDir "${params.output}/ASSEMBLY/UNICYCLER/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}", mode: 'copy'

	input:
		tuple path(pearPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		tuple path("*"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.assembly && !params.skip_assembly && params.run_unicycler) || params.run_all

	script:
		if (reads2 == "-")
			"""
			unicycler -s ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.SE.fq.gz -o . -t $task.cpus --no_correct --mode ${params.assembly_unicycler_mode}
			mv assembly.fasta ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.assembly.fasta
			"""
		else
			"""
			unicycler -1 ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE1.fq.gz -2 ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MEPE2.fq.gz -s ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}/${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.MESE.fq.gz -o . -t $task.cpus --no_correct --mode ${params.assembly_unicycler_mode}
			mv assembly.fasta ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.assembly.fasta
			"""
}

process PROKKA {
	conda baseDir + '/env/snpless-assembly-prokka.yml'
	tag "PROKKA on ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"
	cpus params.prokka_threads

	publishDir "${params.output}/ASSEMBLY/PROKKA", mode: 'copy'

	input:
		tuple path(unicyclerPath), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(proteins)

	output:
		tuple path("${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}"), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	when:
		(params.assembly && !params.skip_assembly && params.run_prokka) || params.run_all

	script:
		"""
		prokka --proteins ${proteins} --outdir ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} --prefix ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType} ${sampleId}_${sampleReplicate}_${sampleTimepoint}_${sampleType}.assembly.fasta
		"""
}
