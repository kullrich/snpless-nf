process SNPEFFCREATEDB {
	conda baseDir + '/env/snpless-annotation-snpeff.yml'
	tag "SNPEFFCREATEDB on ${reference}"
	cpus params.snpeff_threads

	publishDir "${params.output}/ANNOTATION/REFERENCE", mode: 'symlink'

	input:
		path(snpeffDir)
		path(reference)
		path(gff3)

	when:
		(params.annotation && params.run_snpeff && !params.skip_snpeff) || params.run_all

	script:
	    if( params.annotation_snpeff_type == 'gff3' )
		    """
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.genome : ${params.annotation_reference_name}" >> ${snpeffDir}/snpEff.config
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.codonTable : ${params.annotation_reference_codon_table}" >> ${snpeffDir}/snpEff.config
		    mkdir -p ${snpeffDir}/data
		    mkdir -p ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}
		    cp ${reference} ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/sequences.fa
		    cp ${gff3} ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    gzip -f ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    snpEff build -c ${snpeffDir}/snpEff.config -gff3 -v ${params.annotation_reference_name}${params.annotation_reference_version}
		    """
		else if( params.annotation_snpeff_type == 'gtf' )
		    """
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.genome : ${params.annotation_reference_name}" >> ${snpeffDir}/snpEff.config
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.codonTable : ${params.annotation_reference_codon_table}" >> ${snpeffDir}/snpEff.config
		    mkdir -p ${snpeffDir}/data
		    mkdir -p ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}
		    cp ${reference} ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/sequences.fa
		    cp ${gff3} ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    gzip -f ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gff
		    snpEff build -c ${snpeffDir}/snpEff.config -gtf22 -v ${params.annotation_reference_name}${params.annotation_reference_version}
		    """
		else if( params.annotation_snpeff_type == 'genbank' )
		    """
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.genome : ${params.annotation_reference_name}" >> ${snpeffDir}/snpEff.config
		    echo "${params.annotation_reference_name}${params.annotation_reference_version}.codonTable : ${params.annotation_reference_codon_table}" >> ${snpeffDir}/snpEff.config
		    mkdir -p ${snpeffDir}/data
		    mkdir -p ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}
		    cp ${proteins} ${snpeffDir}/data/${params.annotation_reference_name}${params.annotation_reference_version}/genes.gbk
		    snpEff build -c ${snpeffDir}/snpEff.config -genbank -v ${params.annotation_reference_name}${params.annotation_reference_version}
		    """
		else
		    error "Invalid annotation type: ${params.annotation_snpeff_type}"
}

