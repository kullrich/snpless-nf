process FREEBAYESBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/FREEBAYES/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.freebayes.*.vcf", emit: freebayes_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_freebayes && params.run_breseq && !params.skip_breseq && !params.skip_freebayes) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		freebayes -f ${mappingPath}/reference.fasta ${params.snpcalling_freebayes_options} -L bamlist.txt | vcffilter ${params.snpcalling_freebayes_filter_options} > BRESEQ.freebayes.vcffilter.vcf
		vt normalize -r ${mappingPath}/reference.fasta -o BRESEQ.freebayes.normalized.vcf BRESEQ.freebayes.vcffilter.vcf
		vt decompose -o BRESEQ.freebayes.decomposed.vcf BRESEQ.freebayes.normalized.vcf
		"""
}

process FREEBAYESMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/FREEBAYES/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.freebayes.*.vcf", emit: freebayes_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_freebayes && params.run_minimap2 && !params.skip_minimap2 && !params.skip_freebayes) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		freebayes -f ${reference} ${params.snpcalling_freebayes_options} -L bamlist.txt | vcffilter ${params.snpcalling_freebayes_filter_options} > MINIMAP2.freebayes.vcffilter.vcf
		vt normalize -r ${reference} -o MINIMAP2.freebayes.normalized.vcf MINIMAP2.freebayes.vcffilter.vcf
		vt decompose -o MINIMAP2.freebayes.decomposed.vcf MINIMAP2.freebayes.normalized.vcf
		"""
}

process FREEBAYESBWA {
	conda baseDir + '/env/snpless-snpcalling-freebayes.yml'
	tag "FREEBAYES on ${mappingPath}"
	cpus params.freebayes_threads

	publishDir "${params.output}/SNPCALLING/FREEBAYES/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.freebayes.*.vcf", emit: freebayes_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_freebayes && params.run_bwa && !params.skip_bwa && !params.skip_freebayes) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		freebayes -f ${reference} ${params.snpcalling_freebayes_options} -L bamlist.txt | vcffilter ${params.snpcalling_freebayes_filter_options} > BWA.freebayes.vcffilter.vcf
		vt normalize -r ${reference} -o BWA.freebayes.normalized.vcf BWA.freebayes.vcffilter.vcf
		vt decompose -o BWA.freebayes.decomposed.vcf BWA.freebayes.normalized.vcf
		"""
}

process BCFTOOLSBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-bcftools.yml'
	tag "BCFTOOLS on ${mappingPath}"
	cpus params.bcftools_threads

	publishDir "${params.output}/SNPCALLING/BCFTOOLS/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.bcftools.*.vcf", emit: bcftools_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_bcftools && params.run_breseq && !params.skip_breseq && !params.skip_bcftools) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${mappingPath}/reference.fasta -b bamlist.txt | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > BRESEQ.bcftools.varfilter.vcf
		vt normalize -r ${mappingPath}/reference.fasta -o BRESEQ.bcftools.normalized.vcf BRESEQ.bcftools.varfilter.vcf
		vt decompose -o BRESEQ.bcftools.decomposed.vcf BRESEQ.bcftools.normalized.vcf
		"""
}

process BCFTOOLSMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-bcftools.yml'
	tag "BCFTOOLS on ${mappingPath}"
	cpus params.bcftools_threads

	publishDir "${params.output}/SNPCALLING/BCFTOOLS/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.bcftools.*.vcf", emit: bcftools_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_bcftools && params.run_minimap2 && !params.skip_minimap2 && !params.skip_bcftools) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${reference} -b bamlist.txt | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > MINIMAP2.bcftools.varfilter.vcf
		vt normalize -r ${reference} -o MINIMAP2.bcftools.normalized.vcf MINIMAP2.bcftools.varfilter.vcf
		vt decompose -o MINIMAP2.bcftools.decompose.vcf MINIMAP2.bcftools.varfilter.vcf
		"""
}

process BCFTOOLSBWA {
	conda baseDir + '/env/snpless-snpcalling-bcftools.yml'
	tag "BCFTOOLS on ${mappingPath}"
	cpus params.bcftools_threads

	publishDir "${params.output}/SNPCALLING/BCFTOOLS/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.bcftools.*.vcf", emit: bcftools_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_bcftools && params.run_bwa && !params.skip_bwa && !params.skip_bcftools) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		bcftools mpileup ${params.snpcalling_bcftools_mpileup_options} -f ${reference} -b bamlist.txt | bcftools call ${params.snpcalling_bcftools_call_options} | vcfutils.pl varFilter ${params.snpcalling_bcftools_varfilter_options} > BWA.bcftools.varfilter.vcf
		vt normalize -r ${reference} -o BWA.bcftools.normalized.vcf BWA.bcftools.varfilter.vcf
		vt decompose -o BWA.bcftools.decompose.vcf BWA.bcftools.normalized.vcf
		"""
}

process LOFREQBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-lofreq.yml'
	tag "LOFREQBRESEQ on ${sampleLong}"
	cpus params.lofreq_threads

	publishDir "${params.output}/SNPCALLING/LOFREQ/BRESEQ", mode: 'symlink'

	input:
		path(postbreseqDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)

	output:
		path("${sampleLong}"), emit: lofreq_breseqDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: lofreq_breseqSamples
		path "${sampleLong}/${sampleLong}.breseq.lofreq.vcf", emit: lofreq_vcf

	when:
		(params.snpcalling && params.run_lofreq && params.run_breseq && !params.skip_breseq && !params.skip_lofreq) || params.run_all

	script:
		"""
		lofreq faidx ${sampleLong}/data/reference.fasta
		lofreq indelqual --dindel -f ${sampleLong}/data/reference.fasta ${sampleLong}/data/${sampleLong}.breseq.sorted.bam > ${sampleLong}/data/${sampleLong}.breseq.sorted.indel.bam
		lofreq index ${sampleLong}/data/${sampleLong}.breseq.sorted.indel.bam
		lofreq call-parallel --pp-threads $task.cpus ${params.snpcalling_lofreq_options} -f ${sampleLong}/data/reference.fasta -o ${sampleLong}/${sampleLong}.breseq.lofreq.vcf ${sampleLong}/data/${sampleLong}.breseq.sorted.indel.bam
		"""
}

process LOFREQMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-lofreq.yml'
	tag "LOFREQMINIMAP2 on ${sampleLong}"
	cpus params.lofreq_threads

	publishDir "${params.output}/SNPCALLING/LOFREQ/MINIMAP2", mode: 'symlink'

	input:
		path(postminimap2Dir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)
		path(reference)

	output:
		path("${sampleLong}"), emit: lofreq_minimap2Dir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: lofreq_minimap2Samples
		path "${sampleLong}/${sampleLong}.minimap2.lofreq.vcf", emit: lofreq_vcf

	when:
		(params.snpcalling && params.run_lofreq && params.run_minimap2 && !params.skip_minimap2 && !params.skip_lofreq) || params.run_all

	script:
		"""
		lofreq faidx ${reference}
		lofreq indelqual --dindel -f ${reference} ${sampleLong}/${sampleLong}.minimap2.sorted.bam > ${sampleLong}/${sampleLong}.minimap2.sorted.indel.bam
		lofreq index ${sampleLong}/${sampleLong}.minimap2.sorted.indel.bam
		lofreq call-parallel --pp-threads $task.cpus ${params.snpcalling_lofreq_options} -f ${reference} -o ${sampleLong}/${sampleLong}.minimap2.lofreq.vcf ${sampleLong}/${sampleLong}.minimap2.sorted.indel.bam
		"""
}

process LOFREQBWA {
	conda baseDir + '/env/snpless-snpcalling-lofreq.yml'
	tag "LOFREQBWA on ${sampleLong}"
	cpus params.lofreq_threads

	publishDir "${params.output}/SNPCALLING/LOFREQ/BWA", mode: 'symlink'

	input:
		path(postbwaDir)
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2)
		path(reference)

	output:
		path("${sampleLong}"), emit: lofreq_bwaDir
		tuple val(sampleLong), val(sampleId), val(sampleReplicate), val(sampleTimepoint), val(sampleType), val(reads1), val(reads2), emit: lofreq_bwaSamples
		path "${sampleLong}/${sampleLong}.bwa.lofreq.vcf", emit: lofreq_vcf

	when:
		(params.snpcalling && params.run_lofreq && params.run_bwa && !params.skip_bwa && !params.skip_lofreq) || params.run_all

	script:
		"""
		lofreq faidx ${reference}
		lofreq indelqual --dindel -f ${reference} ${sampleLong}/${sampleLong}.bwa.sorted.bam > ${sampleLong}/${sampleLong}.bwa.sorted.indel.bam
		lofreq index ${sampleLong}/${sampleLong}.bwa.sorted.indel.bam
		lofreq call-parallel --pp-threads $task.cpus ${params.snpcalling_lofreq_options} -f ${reference} -o ${sampleLong}/${sampleLong}.bwa.lofreq.vcf ${sampleLong}/${sampleLong}.bwa.sorted.indel.bam
		"""
}

process VARSCANBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "VARSCANBRESEQ on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.varscan.*.vcf", emit: varscan_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_varscan && params.run_breseq && !params.skip_breseq && !params.skip_varscan) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${mappingPath}/reference.fasta -b bamlist.txt | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > BRESEQ.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${mappingPath}/reference.fasta -b bamlist.txt | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > BRESEQ.varscan.indel.vcf
		vt normalize -r ${mappingPath}/reference.fasta -o BRESEQ.varscan.snp.normalized.vcf BRESEQ.varscan.snp.vcf
		vt decompose -o BRESEQ.varscan.snp.decomposed.vcf BRESEQ.varscan.snp.normalized.vcf
		vt normalize -r ${mappingPath}/reference.fasta -o BRESEQ.varscan.indel.normalized.vcf BRESEQ.varscan.indel.vcf
		vt decompose -o BRESEQ.varscan.indel.decompose.vcf BRESEQ.varscan.indel.normalized.vcf
		"""
}

process VARSCANMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "VARSCANMINIMAP2 on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.varscan.*.vcf", emit: varscan_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_varscan && params.run_minimap2 && !params.skip_minimap2 && !params.skip_varscan) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} -b bamlist.txt | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > MINIMAP2.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} -b bamlist.txt | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > MINIMAP2.varscan.indel.vcf
		vt normalize -r ${reference} -o MINIMAP2.varscan.snp.normalized.vcf MINIMAP2.varscan.snp.vcf
		vt decompose -o MINIMAP2.varscan.snp.decompose.vcf MINIMAP2.varscan.snp.normalized.vcf
		vt normalize -r ${reference} -o MINIMAP2.varscan.indel.normalized.vcf MINIMAP2.varscan.indel.vcf
		vt decompose -o MINIMAP2.varscan.indel.decompose.vcf MINIMAP2.varscan.indel.normalized.vcf
		"""
}

process VARSCANBWA {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "VARSCANBWA on ${mappingPath}"
	cpus params.varscan_threads

	publishDir "${params.output}/SNPCALLING/VARSCAN/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.varscan.*.vcf", emit: varscan_vcf
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_varscan && params.run_bwa && !params.skip_bwa && !params.skip_varscan) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} -b bamlist.txt | varscan mpileup2snp ${params.snpcalling_varscan_snp_options} > BWA.varscan.snp.vcf
		samtools mpileup ${params.snpcalling_varscan_mpileup_options} -f ${reference} -b bamlist.txt | varscan mpileup2indel ${params.snpcalling_varscan_indel_options} > BWA.varscan.indel.vcf
		vt normalize -r ${reference} -o BWA.varscan.snp.normalized.vcf BWA.varscan.snp.vcf
		vt decompose -o BWA.varscan.snp.decompose.vcf BWA.varscan.snp.normalized.vcf
		vt normalize -r ${reference} -o BWA.varscan.indel.normalized.vcf BWA.varscan.indel.vcf
		vt decompose -o BWA.varscan.indel.decompose.vcf BWA.varscan.indel.vcf
		"""
}

process MPILEUPBRESEQ {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "MPILEUPBRESEQ on ${mappingPath}"
	cpus params.mpileup_threads

	publishDir "${params.output}/SNPCALLING/MPILEUP/BRESEQ", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)

	output:
		path "BRESEQ.mpileup.txt", emit: mpileup_txt
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_mpileup && params.run_breseq && !params.skip_breseq && !params.skip_mpileup) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		samtools mpileup ${params.snpcalling_mpileup_options} -f ${mappingPath}/reference.fasta -b bamlist.txt | parse_mpileup.py > BRESEQ.mpileup.txt

		"""
}

process MPILEUPMINIMAP2 {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "MPILEUPMINIMAP2 on ${mappingPath}"
	cpus params.mpileup_threads

	publishDir "${params.output}/SNPCALLING/MPILEUP/MINIMAP2", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "MINIMAP2.mpileup.txt", emit: mpileup_txt
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_mpileup && params.run_minimap2 && !params.skip_minimap2 && !params.skip_mpileup) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		samtools mpileup ${params.snpcalling_mpileup_options} -f ${reference} -b bamlist.txt | parse_mpileup.py > MINIMAP2.mpileup.txt

		"""
}

process MPILEUPBWA {
	conda baseDir + '/env/snpless-snpcalling-varscan.yml'
	tag "MPILEUPBWA on ${mappingPath}"
	cpus params.mpileup_threads

	publishDir "${params.output}/SNPCALLING/MPILEUP/BWA", mode: 'symlink'

	input:
		val(bams)
		path(mappingPath)
		path(reference)

	output:
		path "BWA.mpileup.txt", emit: mpileup_txt
		path "bamlist.txt", emit: bamlist

	when:
		(params.snpcalling && params.run_mpileup && params.run_bwa && !params.skip_bwa && !params.skip_mpileup) || params.run_all

	script:
		"""
		sortfiles.py ${mappingPath} bam _ > bamlist.txt
		samtools mpileup ${params.snpcalling_mpileup_options} -f ${reference} -b bamlist.txt | parse_mpileup.py > BWA.mpileup.txt

		"""
}

process GDCOMPARE {
	conda baseDir + '/env/snpless-mapping-breseq.yml'
	tag "GDCOMPARE on ${mappingPath}"
	cpus params.gdtools_threads

	publishDir "${params.output}/SNPCALLING/GDCOMPARE/", mode: 'symlink'

	input:
		val(gds)
		path(mappingPath)
		path(reference)
		path(gff3)

	output:
		path "BRESEQ.annotate.*", emit: gdcompare

	when:
		(params.snpcalling && params.run_gdcompare && params.run_breseq && !params.skip_breseq && !params.skip_gdcompare) || params.run_all

	script:
		"""
		(cat ${gff3}; echo "##FASTA"; cat ${reference}) > reference.gff3
		gdtools ANNOTATE -o BRESEQ.annotate.html -r reference.gff3 ${mappingPath}/*.gd
		gdtools ANNOTATE -f TSV -o BRESEQ.annotate.tsv -r reference.gff3 ${mappingPath}/*.gd
		"""
}
