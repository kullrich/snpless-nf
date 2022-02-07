process SNPEFFCREATEDB {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFCREATEDB on ${reference}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/REFERENCE", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(reference)
		path(gff3)

	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
	    if( params.annotation_snpeff_type == 'gff3' )
		    """
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.genome : ${params.annotation_reference_name}" >> ${snpeffOutDir}/snpEff.config
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.codonTable : ${params.annotation_reference_codon_table}" >> ${snpeffOutDir}/snpEff.config
		    mkdir -p ${snpeffOutDir}/data
		    mkdir -p ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}
		    cp ${reference} ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/sequences.fa
		    cp ${gff3} ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    gzip -f ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    snpEff build -c ${snpeffOutDir}/snpEff.config -gff3 -v ${params.annotation_reference_name}${params.annotation_reference_version}
		    """
		else if( params.annotation_snpeff_type == 'gtf' )
		    """
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.genome : ${params.annotation_reference_name}" >> ${snpeffOutDir}/snpEff.config
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.codonTable : ${params.annotation_reference_codon_table}" >> ${snpeffOutDir}/snpEff.config
		    mkdir -p ${snpeffOutDir}/data
		    mkdir -p ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}
		    cp ${reference} ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/sequences.fa
		    cp ${gff3} ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    gzip -f ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    snpEff build -c ${snpeffOutDir}/snpEff.config -gtf22 -v ${params.annotation_reference_name}${params.annotation_reference_version}
		    """
		else if( params.annotation_snpeff_type == 'genbank' )
		    """
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.genome : ${params.annotation_reference_name}" >> ${snpeffOutDir}/snpEff.config
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.codonTable : ${params.annotation_reference_codon_table}" >> ${snpeffOutDir}/snpEff.config
		    mkdir -p ${snpeffOutDir}/data
		    mkdir -p ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}
		    cp ${proteins} ${snpeffOutDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gbk
		    snpEff build -c ${snpeffOutDir}/snpEff.config -genbank -v ${params.annotation_reference_name}${params.annotation_reference_version}
		    """
		else
		    error "Invalid annotation type: ${params.annotation_snpeff_type}"
}


process SNPEFFANNOTATEFREEBAYESBRESEQ {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/FREEBAYES/BRESEQ", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEFREEBAYESMINIMAP2 {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/FREEBAYES/MINIMAP2", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEFREEBAYESBWA {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/FREEBAYES/BWA", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEBCFTOOLSBRESEQ {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/BCFTOOLS/BRESEQ", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEBCFTOOLSMINIMAP2 {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/BCFTOOLS/MINIMAP2", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEBCFTOOLSBWA {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/BCFTOOLS/BWA", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEVARSCANBRESEQ {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/VARSCAN/BRESEQ", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEVARSCANMINIMAP2 {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/VARSCAN/MINIMAP2", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}

process SNPEFFANNOTATEVARSCANBWA {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFANNOTATE on ${vcfFile}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/VARSCAN/BWA", mode: 'symlink'

	input:
		path(snpeffOutDir)
		path(vcfFile)

	output:
		path "*.vcf", emit: annotation_vcf


	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
		"""
		snpEff -s ${vcfFile}.snpeff.html -csvStats ${vcfFile}.snpeff.csv -c ${snpeffOutDir}/snpEff.config -v ${params.annotation_reference_name}${params.annotation_reference_version} ${vcfFile} > ${vcfFile}.snpeff.vcf
		"""
}
