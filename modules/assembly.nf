process UNICYCLER {
	conda baseDir + '/env/snpless-assembly-unicycler.yml'
	tag "UNICYCLER on ${sampleLong}"
	cpus params.unicycler_threads

	publishDir "${params.output}/ASSEMBLY/UNICYCLER/${sampleLong}", mode: 'copy'

	input:
		path(pearDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		path("*"), emit: unicyclerFiles
		path("${sampleLong}"), emit: unicyclerDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: unicyclerSamples

	when:
		(params.assembly && !params.skip_assembly && params.run_unicycler) || params.run_all

	script:
		if (reads2 == "-")
			"""
			unicycler -s ${sampleLong}/${sampleLong}.SE.fq.gz -o . -t $task.cpus --no_correct --mode ${params.assembly_unicycler_mode}
			mv assembly.fasta ${sampleLong}.assembly.fasta
			"""
		else
			"""
			unicycler -1 ${sampleLong}/${sampleLong}.MEPE1.fq.gz -2 ${sampleLong}/${sampleLong}.MEPE2.fq.gz -s ${sampleLong}/${sampleLong}.MESE.fq.gz -o . -t $task.cpus --no_correct --mode ${params.assembly_unicycler_mode}
			mv assembly.fasta ${sampleLong}.assembly.fasta
			"""
}

process PROKKA {
	conda baseDir + '/env/snpless-assembly-prokka.yml'
	tag "PROKKA on ${sampleLong}"
	cpus params.prokka_threads

	publishDir "${params.output}/ASSEMBLY/PROKKA", mode: 'copy'

	input:
		path(unicyclerFiles)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), path(proteins)

	output:
		path("${sampleLong}"), emit: prokkaDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: prokkaSamples

	when:
		(params.assembly && !params.skip_assembly && params.run_prokka) || params.run_all

	script:
		"""
		prokka --proteins ${proteins} --outdir ${sampleLong} --prefix ${sampleLong} ${sampleLong}.assembly.fasta
		"""
}
