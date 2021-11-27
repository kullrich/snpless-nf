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
        --input                         input table <tab> separated (sampleId;sampleReplicate;sampleTimepoint;reads1,reads2)
        --output                        outdir
        --reference                     path to a fasta file containing the reference
        --gff3                          path to a gff3 file containing the reference annotations
        --proteins                      path to a gbff file containing the reference annotations

    Options: 
        --run_all                       run all: QC, GENMAP, ASSEMBLY, MAPPING, SNPCALLING, SVCALLING, MERGING, ANNOTATION

    Options: QC
        --qc                            run QC

    Options: QC - FASTQC
        --run_fastqc                    run FASTQC
        --fastqc_threads                4
        --qc_fastqc_kmers               length of kmer
        --qc_fastqc_nogroup             disable grouping of bases for reads >50bp
        --qc_fastqc_quiet               supress all progress output on stdout and only report errors

    Options: QC - TRIM READS
        --run_trim                      run TRIM
        --trim_threads                  4
        --qc_trim_adapter_file          path to a fasta file containing all the adapters
        --qc_trim_use_adapter           specify if adapter file should be used
        --qc_clip_options               seedMismatches:palindromeClipThreshold:simpleClipThreshold:
        --qc_trim_options               "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:50"
        --qc_trim_quiet                 supress all progress output on stdout and only report errors

    Options: QC - PEAR READS
        --run_pear                      run PEAR
        --pear_threads                  4
        --qc_pear_options               "-p 0.01 -v 10 -m 999999 -n 50 -t 1 -q 0 -u 1 -g1 -s 2"


    Options: GENMAP
        --genmap                        run GENMAP

    Options: GENMAP - INDEX/MAP
        --run_genmap                    run GENMAP INDEX/MAP
        --genmap_threads                8
        --genmap_kmers                  30
        --genmap_errors                 2
        --genmap_outputs                "-t -w -b"

    Options: ASSEMBLY
        --assembly                      run ASSEMBLY

    Options: ASSEMBLY - UNICYCLER
        --run_unicycler                 run UNICYCLER
        --unicycler_threads             8
        --assembly_unicycler_mode       unicycler bridging mode

    Options: ASSEMBLY - PROKKA
        --run_prokka                    run PROKKA
        --prokka_threads                8

    Options: MAPPING
        --mapping                       run MAPPING

    Options: MAPPING - BRESEQ
        --run_breseq                    run BRESEQ
        --breseq_threads                8
        --mapping_breseq_p              The sample is not clonal. Predict polymorphic (mixed) mutations.
        --mapping_breseq_options        "-m 20 -b 15"
        --mapping_breseq_coverage       "--min-BQ 3 --min-MQ 10"

    Options: MAPPING - MINIMAP2
        --run_minimap2                  run MINIMAP2
        --minimap2_threads              8
        --mapping_minimap2_options      "--sam-hit-only --secondary=yes -ax sr"
        --mapping_minimap2_samblaster   remove duplicates with samblaster
        --mapping_minimap2_coverage     "--min-BQ 3 --min-MQ 10"

    Options: MAPPING - BWA
        --run_bwa                       run BWA
        --bwa_threads                   8
        --mapping_bwa_options           "-M"
        --mapping_bwa_samblaster        remove duplicates with samblaster
        --mapping_bwa_coverage          "--min-BQ 3 --min-MQ 10"

"""
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// INCLUDES
include {FASTQC; TRIM; PEAR} from './modules/qc' params(params)
include {GENMAP} from './modules/genmap' params(params)
include {UNICYCLER; PROKKA} from './modules/assembly' params(params)
include {BRESEQ; MINIMAP2; BWA; POSTBRESEQ} from './modules/mapping' params(params)

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
"""
                paramsFile  = file(params.input).readLines().each{ println it }
            }
            // SET SAMPLES FROM INPUT
            channel
                .fromPath(params.input)
                .splitCsv(header:false,sep:'\t')
                .map{row->tuple(row[0],row[1],row[2],row[3],row[4],file(row[3]),file(row[4]),file(params.qc_trim_adapter_file))}
                .set{samples}
            // run all
            // QC
            // PROCESS FASTQC
            FASTQC(samples)
            // FASTQC.out.view()
            // PROCESS TRIM READS
            TRIM(samples)
            // TRIM.out.view()
            // PROCESS PEAR READS
            PEAR(TRIM.out)
            // PEAR.out.view()
            // GENMAP
            // PROCESS GENMAP
            GENMAP(file(params.reference))
            // GENMAP.out.view()
            // ASSEMBLY
            // PROCESS UNICYCLER
            UNICYCLER(PEAR.out)
            // UNICYCLER.out.view()
            // PROCESS PROKKA
            PROKKA(UNICYCLER.out.map{it + [file(params.proteins)]})
            // PROKKA.out.view()
            // MAPPING
            // PROCESS BRESEQ
            BRESEQ(PEAR.out.map{it + [file(params.proteins)]})
            // BRESEQ.out.view()
            POSTBRESEQ(BRESEQ.out)
            // POSTBRESEQ.out.view()
            // PROCESS MINIMAP2
            MINIMAP2(PEAR.out.map{it + [file(params.reference)]})
            // MINIMAP2.out.view()
            // PROCESS BWA
            BWA(PEAR.out.map{it + [file(params.reference)]})
            // BWA.out.view()

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
