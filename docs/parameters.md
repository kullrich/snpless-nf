## Parameters

## Basic options

:question: `--help`

`--version` Print version.

:bug: `--debug` Turns on debug mode.

## :arrow_forward: Input/Output options

Define where the pipeline should find input data and save output data.

:floppy_disk: `--input <ARG>` Path to tab-separated file containing information about the samples in the experiment. :pushpin: `mandatory`

:open_file_folder: `--output <results>` Path to the ouptut directory where the results will be saved. `default: results`

:floppy_disk: `--reference <ARG>` Path to a FASTA reference file. :pushpin: `mandatory`

:floppy_disk: `--gff3 <ARG>` Path to a GFF3 file for the given reference. :pushpin: `mandatory`

:floppy_disk: `--proteins <ARG>` Path to a GBFF (Genomic GenBank format) reference file. :pushpin: `mandatory`

## :information_source: General options

:hourglass_flowing_sand: `--run_all` Turns on all pipeline steps. `default: true` 

:mag: `--qc` Turns on QC step. `default: true`

:earth_americas: `--genmap` Turns on GENMAP step. `default: true`

:wrench: `--assembly` Turns on ASSEMBLY step. `default: true`

:dart: `--mapping` Turns on MAPPING step. `default: true`

:dango: `--snpcalling` Turns on SNPCALLING step. `default: true`

:oden: `--svcalling` Turns on SVCALLING step. `default: true`

## :fast_forward: Skip steps

Skip any of the mentioned steps.

Some of the steps in the pipeline can be executed optionally. If you specify specific steps to be skipped, there won't be any output related to these modules.

:fast_forward: `--skip_genmap` Turns off GENMAP mappability calculation for the reference. This information is normally used to exclude some regions from SVCALLNG.

:fast_forward: `--skip_assembly` Turns off de-novo ASSEMBLY with UNICYCLER and subsequent protein prediction with PROKKA.

:fast_forward: `--skip_breseq` Turns off MAPPPING with BRESEQ and all related downstream analysis of BRESEQ output.

:fast_forward: `--skip_minimap2` Turns off MAPPING with MINIMAP2 and all related downstream analysis of MINIMAP2 output.

:fast_forward: `--skip_bwa` Turns off MAPPING with BWA and all related downstream analysis of BWA output.

:fast_forward: `--skip_freebayes` Turns off SNPCALLING with FREEBAYES and all related downstream analysis of FREEBAYES output.

:fast_forward: `--skip_bcftools` Turns off SNPCALLING with BCFTOOLS and all related downstream analysis of BCFTOOLS output.

:fast_forward: `--skip_lofreq` Turns off SNPCALLING with LOFREQ and all related downstream analysis of LOFREQ output.

:fast_forward: `--skip_varscan` Turns off SNPCALLING with VARSCAN and all related downstream analysis of VARSCAN output.

:fast_forward: `--skip_gdcompare` Turns off SNPCALLING with GDCOMAPRE and all related downstream analysis of GDCOMPARE output.

:fast_forward: `--skip_pindel` Turns off SNPCALLING with PINDEL and all related downstream analysis of PINDEL output.

:fast_forward: `--skip_gridss` Turns off SNPCALLING with GRIDSS and all related downstream analysis of GRIDSS output.

## QC options

Options to adjust read quality.

### QC - FASTQC

`--fastqc_threads` 

`--qc_fastqc_kmers` 

`--qc_fastqc_nogroup`

`--qc_fastqc_quiet`

### QC - TRIM

`--trim_threads`

:scissors: `--qc_trim_adapter_file` path to a fasta file containing all the adapters

:scissors: `--qc_trim_use_adapter` specify if adapter file should be used

:scissors: `--qc_clip_options` seedMismatches:palindromeClipThreshold:simpleClipThreshold:

:scissors: `--qc_trim_options` "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:50"

:scissors: `--qc_trim_quiet` supress all progress output on stdout and only report errors

## GENMAP options

## ASSEMBLY options

## MAPPING options

## SNPCALLING options

## SVCALLING options

