#!/usr/bin/env nextflow

/*
Copyright (c) 2021, Kristian K Ullrich, Loukas Theodosiou.
*/


/*
    ========
    snpless-nf - A Nextflow pipeline for time-course analysis with bacterial NGS whole-genome data.
    ========
    Kristian K Ullrich
    Loukas Theodosiou
    --------
    Github : https://
*/

nextflow.enable.dsl = 2

/*
    ========
    VERSION
    ========
*/
if(params.version){
    println """\
========
SNPLESS
========
~ version ${workflow.manifest.version}
"""
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

/*
    ========
    HELP
    ========
*/
if(params.help){
    println """\
========
SNPLESS
========
~ version ${workflow.manifest.version}

Usage:
    nextflow run snpless-nf [OPTIONS]...

    Options: INPUT/OUTPUT
        --input                                 input table <tab> separated (sampleId;sampleReplicate;sampleTimepoint;sampleType;reads1,reads2)
        --output                                outdir
        --reference                             path to a fasta file containing the reference
        --gff3                                  path to a gff3 file containing the reference annotations
        --gtf                                   path to a gtf file containing the reference annotations
        --genbank                               path to a gbff file containing the reference annotations
        --proteins                              path to a gbff file containing the reference annotations

    Options: 
        --run_all                               run all: QC, GENMAP, ASSEMBLY, MAPPING, SNPCALLING, SVCALLING, MERGING, ANNOTATION

    Options: SKIP STEPS
        --skip_genmap
        --skip_assembly
        --skip_breseq
        --skip_minimap2
        --skip_bwa
        --skip_freebayes
        --skip_bcftools
        --skip_lofreq
        --skip_varscan
        --skip_gdcompare
        --skip_pindel
        --skip_gridss
        --skip_snpeff

    Options: QC
        --qc                                    run QC

    Options: QC - FASTQC
        --run_fastqc                            run FASTQC
        --fastqc_threads                        4
        --qc_fastqc_kmers                       length of kmer
        --qc_fastqc_nogroup                     disable grouping of bases for reads >50bp
        --qc_fastqc_quiet                       supress all progress output on stdout and only report errors

    Options: QC - TRIM READS
        --run_trim                              run TRIM
        --trim_threads                          4
        --qc_trim_adapter_file                  path to a fasta file containing all the adapters
        --qc_trim_use_adapter                   specify if adapter file should be used
        --qc_clip_options                       seedMismatches:palindromeClipThreshold:simpleClipThreshold:
        --qc_trim_options                       "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:50"
        --qc_trim_quiet                         supress all progress output on stdout and only report errors

    Options: QC - PEAR READS
        --run_pear                              run PEAR
        --pear_threads                          4
        --qc_pear_options                       "-p 0.01 -v 10 -m 999999 -n 50 -t 1 -q 0 -u 1 -g1 -s 2"


    Options: GENMAP
        --genmap                                run GENMAP

    Options: GENMAP - INDEX/MAP
        --run_genmap                            run GENMAP INDEX/MAP
        --genmap_threads                        8
        --genmap_kmers                          30
        --genmap_errors                         2
        --genmap_outputs                        "-t -w -b"

    Options: ASSEMBLY
        --assembly                              run ASSEMBLY

    Options: ASSEMBLY - UNICYCLER
        --run_unicycler                         run UNICYCLER
        --unicycler_threads                     8
        --assembly_unicycler_mode               unicycler bridging mode

    Options: ASSEMBLY - PROKKA
        --run_prokka                            run PROKKA
        --prokka_threads                        8

    Options: MAPPING
        --mapping                               run MAPPING

    Options: MAPPING - BRESEQ
        --run_breseq                            run BRESEQ
        --breseq_threads                        8
        --mapping_breseq_p                      sample is not clonal. Predict polymorphic (mixed) mutations.
        --mapping_breseq_options                "-m 30 -b 20"
        --mapping_breseq_coverage               "--min-MQ 30 --min-BQ 20"

    Options: MAPPING - MINIMAP2
        --run_minimap2                          run MINIMAP2
        --minimap2_threads                      8
        --mapping_minimap2_options              "--sam-hit-only --secondary=yes -ax sr"
        --mapping_minimap2_samblaster           remove duplicates with samblaster
        --mapping_minimap2_coverage             "--min-MQ 30 --min-BQ 20"

    Options: MAPPING - BWA
        --run_bwa                               run BWA
        --bwa_threads                           8
        --mapping_bwa_options                   "-M"
        --mapping_bwa_samblaster                remove duplicates with samblaster
        --mapping_bwa_coverage                  "--min-MQ 30 --min-BQ 20"

    Options: SNPCALLING
        --snpcalling                            run SNPCALLING

    Options: SNPCALLING - FREEBAYES
        --run_freebayes                         run FREEBAYES
        --freebayes_threads                     1
        --snpcalling_freebayes_options          "--pooled-discrete --min-alternate-fraction 0.05 --min-alternate-count 2 --min-mapping-quality 30 --min-base-quality 20"
        --snpcalling_freebayes_filter_options   "-f 'QUAL > 30'"

    Options: SNPCALLING - BCFFTOOLS
        --run_bcftools                          run BCFTOOLS
        --bcftools_threads                      1
        --snpcalling_bcftools_mpileup_options   "-C 50 -q 30 -Q 20 -F 0.0002 -d 2000 -E -a FORMAT/AD,FORMAT/DP"
        --snpcalling_bcftools_call_options      "-mAv -Ov"
        --snpcalling_bcftools_varfilter_options "-Q 10 -d 5 -D 2000"

    Options: SNPCALLING - GDCOMPARE
        --run_gdcompare                         run GDCOMPARE
        --gdtools_threads                       1

    Options: SNPCALLING - LOFREQ
        --run_lofreq                            run LOFREQ
        --lofreq_threads                        8
        --snpcalling_lofreq_options             "-C 5 -d 2000 -m 30 -q 20 -Q 20 -D --call-indels"

    Options: SNPCALLING - VARSCAN
        --run_varscan                           run VARSCAN
        --varscan_threads                       1
        --snpcalling_varscan_mpileup_options    "-C 50 -q 30 -Q 20 -d 2000 -E"
        --snpcalling_varscan_snp_options        "--min-coverage 5 --min-avg-qual 20 --min-var-freq 0.05 --output-vcf 1"
        --snpcalling_varscan_indel_options      "--min-coverage 5 --min-avg-qual 20 --min-var-freq 0.05 --output-vcf 1"

    Options: SVCALLING                          
        --svcalling                             run SVCALLING

    Options: SVCALLING - PINDEL
        --run_pindel                            run PINDEL
        --pindel_threads                        8
        --svcalling_pindel_sam2pindel_options   "300 tag 0 Illumina-PairEnd"
        --svcalling_pindel_options              "-c ALL"

    Options: SVCALLING - GRIDSS

    Options: ANNOTATION
        --annotation                            run ANNOTATION

    Options: ANNOTATION - SNPEFF
        --run_snpeff                            run SNPEFF
        --annotation_snpeff_type                "gff3"
        --annotation_reference_name             "refname"
        --annotation_reference_version          "0.0"
        --annotation_reference_codon_table      "Bacterial_and_Plant_Plastid"
        --annotation_snpeff_config_file         "data/snpEff.config"

nextflow run snpless-nf --input <samples.tsv> --reference <genome.fna> --gff3 <genome.gff3> --proteins <genome.gbff>

"""
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// INCLUDES
include {FASTQC; TRIM; PEAR} from './modules/qc' params(params)
include {GENMAP} from './modules/genmap' params(params)
include {UNICYCLER; PROKKA} from './modules/assembly' params(params)
include {BRESEQ; MINIMAP2; BWA; POSTBRESEQ; POSTMINIMAP2; POSTBWA} from './modules/mapping' params(params)
include {FREEBAYESBRESEQ; FREEBAYESMINIMAP2; FREEBAYESBWA; BCFTOOLSBRESEQ; BCFTOOLSMINIMAP2; BCFTOOLSBWA; GDCOMPARE; LOFREQBRESEQ; LOFREQMINIMAP2; LOFREQBWA; VARSCANBRESEQ; VARSCANMINIMAP2; VARSCANBWA} from './modules/snpcalling' params(params)
include {PINDELMINIMAP2; PINDELBWA} from './modules/svcalling' params(params)
include {SNPEFFCREATEDB} from './modules/annotation' params(params)

// QC workflow
workflow qc {
    take:
        samples
    main:
        //
        // PROCESS FASTQC
        FASTQC(samples)
        // FASTQC.out.subscribe {println "Got: $it"}
        // FASTQC.out.view()
        // FASTQC.out.subscribe onComplete: {println "FastQC - Done"}
        //
        // PROCESS TRIM READS
        TRIM(samples)
        // TRIM.out.subscribe {println "Got: $it"}
        // TRIM.out.view()
        // TRIM.out.subscribe onComplete: {println "Trimmomatic - Done"}
        //
        // PROCESS PEAR READS
        PEAR(TRIM.out)
        // PEAR.out.subscribe {println "Got: $it"}
        // PEAR.out.view()
        // PEAR.out.subscribe onComplete: {println "Pear - Done"}
    emit:
        fastqc = FASTQC.out
        trim = TRIM.out
        pear = PEAR.out
}

// GENMAP workflow
workflow genmap {
    main:
        // PROCESS GENMAP
        GENMAP(file(params.reference))
        // GENMAP.out.subscribe {println "Got: $it"}
        // GENMAP.out.view()
        // GENMAP.out.subscribe onComplete: {println "GenMap - Done"}
    emit:
        genmap = GENMAP.out
}

// ASSEMBLY workflow
workflow assembly {
    take:
        pear
    main:
        // PROCESS UNICYCLER
        UNICYCLER(pear)
        // UNICYCLER.out.subscribe {println "Got: $it"}
        // UNICYCLER.out.view()
        // UNICYCLER.out.subscribe onComplete: {println "Unicycler - Done"}
        //
        // PROCESS PROKKA
        PROKKA(UNICYCLER.out.map{it + [file(params.proteins)]})
        // PROKKA.out.subscribe {println "Got: $it"}
        // PROKKA.out.view()
        // PROKKA.out.subscribe onComplete: {println "Prokka - Done"}
    emit:
        unicycler = UNICYCLER.out
        prokka = PROKKA.out
}

workflow mapping {
    take:
        pear
        breseqDir
        minimap2Dir
        bwaDir
    main:
        // PROCESS BRESEQ
        BRESEQ(pear.map{it + [file(params.proteins)]})
        // BRESEQ.out.subscribe {println "Got: $it"}
        // BRESEQ.out.view()
        // BRESEQ.out.subscribe onComplete: {println "breseq - Done"}
        POSTBRESEQ(BRESEQ.out)
        // POSTBRESEQ.out.subscribe {println "Got: $it"}
        // POSTBRESEQ.out.view()
        // POSTBRESEQ.out.subscribe onComplete: {println "Postprocess breseq - Done"}
        POSTBRESEQ.out.breseq_mean_coverage.subscribe {it.copyTo(breseqDir)}
        POSTBRESEQ.out.breseq_bam.subscribe {it.copyTo(breseqDir)}
        POSTBRESEQ.out.breseq_bam_index.subscribe {it.copyTo(breseqDir)}
        POSTBRESEQ.out.breseq_bam_reference.subscribe {it.copyTo(breseqDir)}
        POSTBRESEQ.out.breseq_bam_reference_gff3.subscribe {it.copyTo(breseqDir)}
        POSTBRESEQ.out.breseq_vcf.subscribe {it.copyTo(breseqDir)}
        POSTBRESEQ.out.breseq_gd.subscribe {it.copyTo(breseqDir)}
        postbreseq = POSTBRESEQ.out.postbreseq
        breseq_mean_coverage = POSTBRESEQ.out.breseq_mean_coverage.collect()
        breseq_bam = POSTBRESEQ.out.breseq_bam.collect()
        breseq_bam_index = POSTBRESEQ.out.breseq_bam_index.collect()
        breseq_bam_reference = POSTBRESEQ.out.breseq_bam_reference.collect()
        breseq_bam_reference_gff3 = POSTBRESEQ.out.breseq_bam_reference_gff3.collect()
        breseq_vcf = POSTBRESEQ.out.breseq_vcf.collect()
        breseq_gd = POSTBRESEQ.out.breseq_gd.collect()
        //
        // PROCESS MINIMAP2
        MINIMAP2(pear.map{it + [file(params.reference)]})
        // MINIMAP2(PEAR.out.map{it + [file(params.reference)]})
        // MINIMAP2.out.subscribe {println "Got: $it"}
        // MINIMAP2.out.view()
        // MINIMAP2.out.subscribe onComplete: {println "minimap2 - Done"}
        POSTMINIMAP2(MINIMAP2.out)
        // POSTMINIMAP2.out.subscribe {println "Got: $it"}
        // POSTMINIMAP2.out.view()
        // POSTMINIMAP2.out.subscribe onComplete: {println "Postprocess minimap2 - Done"}
        POSTMINIMAP2.out.minimap2_mean_coverage.subscribe {it.copyTo(minimap2Dir)}
        POSTMINIMAP2.out.minimap2_bam.subscribe {it.copyTo(minimap2Dir)}
        POSTMINIMAP2.out.minimap2_bam_index.subscribe {it.copyTo(minimap2Dir)}
        postminimap2 = POSTMINIMAP2.out.postminimap2
        minimap2_mean_coverage = POSTMINIMAP2.out.minimap2_mean_coverage.collect()
        minimap2_bam = POSTMINIMAP2.out.minimap2_bam.collect()
        minimap2_bam_index = POSTMINIMAP2.out.minimap2_bam_index.collect()
        //
        // PROCESS BWA
        BWA(pear.map{it + [file(params.reference)]})
        // BWA.out.subscribe {println "Got: $it"}
        // BWA.out.view()
        // BWA.out.subscribe onComplete: {println "bwa - Done"}
        POSTBWA(BWA.out)
        // POSTBWA.out.subscribe {println "Got: $it"}
        // POSTBWA.out.view()
        // POSTBWA.out.subscribe onComplete: {println "Postprocess bwa - Done"}
        POSTBWA.out.bwa_mean_coverage.subscribe {it.copyTo(bwaDir)}
        POSTBWA.out.bwa_bam.subscribe {it.copyTo(bwaDir)}
        POSTBWA.out.bwa_bam_index.subscribe {it.copyTo(bwaDir)}
        postbwa = POSTBWA.out.postbwa
        bwa_mean_coverage = POSTBWA.out.bwa_mean_coverage.collect()
        bwa_bam = POSTBWA.out.bwa_bam.collect()
        bwa_bam_index = POSTBWA.out.bwa_bam_index.collect()
    emit:
        postbreseq
        breseq_mean_coverage
        breseq_bam
        breseq_bam_index
        breseq_bam_reference
        breseq_bam_reference_gff3
        breseq_vcf
        breseq_gd
        postminimap2
        minimap2_mean_coverage
        minimap2_bam
        minimap2_bam_index
        postbwa
        bwa_mean_coverage
        bwa_bam
        bwa_bam_index
}

workflow snpcalling {
    take:
        breseqDir
        minimap2Dir
        bwaDir
        lofreq_breseqDir
        lofreq_minimap2Dir
        lofreq_bwaDir
        postbreseq
        breseq_mean_coverage
        breseq_bam
        breseq_bam_index
        breseq_bam_reference
        breseq_bam_reference_gff3
        breseq_vcf
        breseq_gd
        postminimap2
        minimap2_mean_coverage
        minimap2_bam
        minimap2_bam_index
        postbwa
        bwa_mean_coverage
        bwa_bam
        bwa_bam_index
    main:
        // PROCESS FREEBAYESBRESEQ
        FREEBAYESBRESEQ(breseq_bam, breseqDir)
        // PROCESS FREEBAYESMINIMAP2
        FREEBAYESMINIMAP2(minimap2_bam, minimap2Dir, file(params.reference))
        // PROCESS FREEBAYESBWA
        FREEBAYESBWA(bwa_bam, bwaDir, file(params.reference))
        // PROCESS BCFTOOLSBRESEQ
        BCFTOOLSBRESEQ(breseq_bam, breseqDir)
        // PROCESS BCFTOOLSMINIMAP2
        BCFTOOLSMINIMAP2(minimap2_bam, minimap2Dir, file(params.reference))
        // PROCESS BCFTOOLSBWA
        BCFTOOLSBWA(bwa_bam, bwaDir, file(params.reference))
        // PROCESS LOFREQBRESEQ
        LOFREQBRESEQ(postbreseq)
        LOFREQBRESEQ.out.lofreq_vcf.subscribe {it.copyTo(lofreq_breseqDir)}
        // PROCESS LOFREQMINIMAP2
        LOFREQMINIMAP2(postminimap2, file(params.reference))
        LOFREQMINIMAP2.out.lofreq_vcf.subscribe {it.copyTo(lofreq_minimap2Dir)}
        // PROCESS LOFREQBWA
        LOFREQBWA(postbwa, file(params.reference))
        LOFREQBWA.out.lofreq_vcf.subscribe {it.copyTo(lofreq_bwaDir)}
        // PROCESS VARSCANBRESEQ
        VARSCANBRESEQ(breseq_bam, breseqDir)
        // PROCESS VARSCANMINIMAP2
        VARSCANMINIMAP2(minimap2_bam, minimap2Dir, file(params.reference))
        // PROCESS VARSCANBWA
        VARSCANBWA(bwa_bam, bwaDir, file(params.reference))
        // PROCESS GDCOMPARE
        GDCOMPARE(breseq_gd, breseqDir, file(params.proteins))
        freebayes_breseq_vcf = FREEBAYESBRESEQ.out.freebayes_vcf
        freebayes_minimap2_vcf = FREEBAYESMINIMAP2.out.freebayes_vcf
        freebayes_bwa_vcf = FREEBAYESBWA.out.freebayes_vcf
        bcftools_breseq_vcf = BCFTOOLSBRESEQ.out.bcftools_vcf
        bcftools_minimap2_vcf = BCFTOOLSMINIMAP2.out.bcftools_vcf
        bcftools_bwa_vcf = BCFTOOLSBWA.out.bcftools_vcf
        varscan_breseq_vcf = VARSCANBRESEQ.out.varscan_vcf
        varscan_minimap2_vcf = VARSCANMINIMAP2.out.varscan_vcf
        varscan_bwa_vcf = VARSCANBWA.out.varscan_vcf
        lofreq_breseq_vcf = LOFREQBRESEQ.out.lofreq_vcf
        lofreq_minimap2_vcf = LOFREQMINIMAP2.out.lofreq_vcf
        lofreq_bwa_vcf = LOFREQBWA.out.lofreq_vcf
        gdcompare = GDCOMPARE.out.gdcompare
    emit:
        freebayes_breseq_vcf
        freebayes_minimap2_vcf
        freebayes_bwa_vcf
        bcftools_breseq_vcf
        bcftools_minimap2_vcf
        bcftools_bwa_vcf
        varscan_breseq_vcf
        varscan_minimap2_vcf
        varscan_bwa_vcf
        lofreq_breseq_vcf
        lofreq_minimap2_vcf
        lofreq_bwa_vcf
        gdcompare
}

workflow svcalling {
    take:
        breseqDir
        minimap2Dir
        bwaDir
        pindel_breseqDir
        pindel_minimap2Dir
        pindel_bwaDir
        postbreseq
        breseq_mean_coverage
        breseq_bam
        breseq_bam_index
        breseq_bam_reference
        breseq_bam_reference_gff3
        breseq_vcf
        breseq_gd
        postminimap2
        minimap2_mean_coverage
        minimap2_bam
        minimap2_bam_index
        postbwa
        bwa_mean_coverage
        bwa_bam
        bwa_bam_index
    main:
        // PROCESS PINDELBRESEQ
        // PINDELBRESEQ(postbreseq)
        // PINDELBRESEQ.out.pindel_files.subscribe {it.each {it.copyTo(pindel_breseqDir)} }
        // PROCESS PINDELMINIMAP2
        PINDELMINIMAP2(postminimap2, file(params.reference))
        PINDELMINIMAP2.out.pindel_files.subscribe {it.each {it.copyTo(pindel_minimap2Dir)} }
        // PROCESS PINDELBWA
        PINDELBWA(postbwa, file(params.reference))
        PINDELBWA.out.pindel_files.subscribe {it.each {it.copyTo(pindel_bwaDir)} }
        // PROCESS GRIDSSBRESEQ
        // GRIDSSBRESEQ(breseq_bam, breseqDir)
        // PROCESS GRIDSSMINIMAP2
        // GRIDSSMINIMAP2(minimap2_bam, minimap2Dir, file(params.reference))
        // PROCESS GRIDSSBWA
        // GRIDSSBWA(bwa_bam, bwaDir, file(params.reference))
}

workflow annotation {
    take:
        snpeffDir
    main:
        // CREATE REFERENCE DB
        file(params.annotation_snpeff_config_file).copyTo(snpeffDir)
        SNPEFFCREATEDB(snpeffDir, file(params.reference), file(params.gff3))
}

// MAIN workflow
workflow{
    main:
        if (params.input){
            if (params.debug){
            // DEBUG - params file
                println """\
========
SNPLESS
========
~ version ${workflow.manifest.version}
========
DEBUG
--------
INPUT: ${params.input}
OUTPUT: ${params.output}
"""
                paramsFile  = file(params.input).readLines().each{ println it }
            }
            // SET SAMPLES FROM INPUT
            channel
                .fromPath(params.input)
                .splitCsv(header:false,sep:'\t')
                .map{row->tuple(row[0],row[1],row[2],row[3],row[4],row[5],file(row[4]),file(row[5]),file(params.qc_trim_adapter_file))}
                .set{samples}
            // run all
            //
            // QC
            qc(samples)
            //
            // GENMAP
            //
            genmap()
            //
            // ASSEMBLY
            //
            assembly(qc.out.pear)
            //
            // MAPPING
            //
            breseqDir = file(params.output+"/MAPPING/BRESEQOUT")
            breseqDir.mkdirs()
            // println breseqDir
            minimap2Dir = file(params.output+"/MAPPING/MINIMAP2OUT")
            minimap2Dir.mkdirs()
            // println minimap2Dir
            bwaDir = file(params.output+"/MAPPING/BWAOUT")
            bwaDir.mkdirs()
            // println bwaDir
            mapping(qc.out.pear, breseqDir, minimap2Dir, bwaDir)
            //
            // SNPCALLING
            //
            lofreq_breseqDir = file(params.output+"/SNPCALLING/LOFREQ/BRESEQ")
            lofreq_breseqDir.mkdirs()
            lofreq_minimap2Dir = file(params.output+"/SNPCALLING/LOFREQ/MINIMAP2")
            lofreq_minimap2Dir.mkdirs()
            lofreq_bwaDir = file(params.output+"/SNPCALLING/LOFREQ/BWA")
            lofreq_bwaDir.mkdirs()
            //
            snpcalling(breseqDir, minimap2Dir, bwaDir, lofreq_breseqDir, lofreq_minimap2Dir, lofreq_bwaDir, mapping.out.postbreseq, mapping.out.breseq_mean_coverage, mapping.out.breseq_bam, mapping.out.breseq_bam_index, mapping.out.breseq_bam_reference, mapping.out.breseq_bam_reference_gff3, mapping.out.breseq_vcf, mapping.out.breseq_gd, mapping.out.postminimap2, mapping.out.minimap2_mean_coverage, mapping.out.minimap2_bam, mapping.out.minimap2_bam_index, mapping.out.postbwa, mapping.out.bwa_mean_coverage, mapping.out.bwa_bam, mapping.out.bwa_bam_index)
            //
            // SVCALLING
            //
            pindel_breseqDir = file(params.output+"/SVCALLING/PINDEL/BRESEQ")
            pindel_breseqDir.mkdirs()
            pindel_minimap2Dir = file(params.output+"/SVCALLING/PINDEL/MINIMAP2")
            pindel_minimap2Dir.mkdirs()
            pindel_bwaDir = file(params.output+"/SVCALLING/PINDEL/BWA")
            pindel_bwaDir.mkdirs()
            //
            svcalling(breseqDir, minimap2Dir, bwaDir, pindel_breseqDir, pindel_minimap2Dir, pindel_bwaDir, mapping.out.postbreseq, mapping.out.breseq_mean_coverage, mapping.out.breseq_bam, mapping.out.breseq_bam_index, mapping.out.breseq_bam_reference, mapping.out.breseq_bam_reference_gff3, mapping.out.breseq_vcf, mapping.out.breseq_gd, mapping.out.postminimap2, mapping.out.minimap2_mean_coverage, mapping.out.minimap2_bam, mapping.out.minimap2_bam_index, mapping.out.postbwa, mapping.out.bwa_mean_coverage, mapping.out.bwa_bam, mapping.out.bwa_bam_index)
            //
            // MERGING
            //

            //
            // ANNOTATION
            //
            snpeffDir = file(params.output+"/ANNOTATION/REFERENCE")
            snpeffDir.mkdirs()
            annotation(snpeffDir)
            //
        }
}

// WORKFLOW TRACING
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${workflow.success && !params.debug ? "(cleared)" : ""}"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""
}
