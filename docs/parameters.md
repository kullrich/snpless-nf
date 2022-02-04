## Parameters

## Basic options

:question: `--help`

`--version` Print version.

:bug: `--debug` Turns on debug mode.

## :arrow_forward: Input/Output options

Define where the pipeline should find input data and save output data.

:floppy_disk: `--input <ARG>` Path to tab-separated file containing information about the samples in the experiment (sampleId;sampleReplicate;sampleTimepoint;sampleType;reads1,reads2). :pushpin: `mandatory`

:open_file_folder: `--output <results>` Path to the ouptut directory where the results will be saved. `default: results`

:floppy_disk: `--reference <ARG>` Path to a FASTA reference file. :pushpin: `mandatory`

:floppy_disk: `--gff3 <ARG>` Path to a GFF3 file for the reference annotations. :pushpin: `mandatory`

:floppy_disk: `--gtf <ARG>` Path to a GTF file for the reference annotations. :pushpin: `mandatory`

:floppy_disk: `--genbank <ARG>` Path to a GBFF (Genomic GenBank format) reference file. :pushpin: `mandatory`

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

## :mag: QC options

Options to adjust read quality.

### QC - FASTQC

`--fastqc_threads` 

:mag: `--qc_fastqc_kmers` 

:mag: `--qc_fastqc_nogroup`

:mag: `--qc_fastqc_quiet`

### QC - TRIM

see here for detailed possible options [http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

`--trim_threads`

:scissors: `--qc_trim_adapter_file` path to a fasta file containing all the adapters

:scissors: `--qc_trim_use_adapter` specify if adapter file should be used

:scissors: `--qc_clip_options` seedMismatches:palindromeClipThreshold:simpleClipThreshold:

:scissors: `--qc_trim_options` "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:50"

:scissors: `--qc_trim_quiet` supress all progress output on stdout and only report errors

### QC - PEAR

see here for detailed possible options [https://cme.h-its.org/exelixis/web/software/pear/doc.html](https://cme.h-its.org/exelixis/web/software/pear/doc.html)

`--pear_threads`

`--qc_pear_options` "-p 0.01 -v 10 -m 999999 -n 50 -t 1 -q 0 -u 1 -g1 -s 2"

## :earth_americas: GENMAP options

`--genmap_threads`

`--genmap_kmers`

`--genmap_errors`

`--genmap_outputs`

## :wrench: ASSEMBLY options

### ASSEMBLY - UNICYCLER

`--unicycler_threads`

`--assembly_unicycler_mode`

### ASSEMBLY - PROKKA

`--prokka_threads`

## :dart: MAPPING options

### MAPPING - BRESEQ

see here for detailed possible options [https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/usage.html#breseq](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/usage.html#breseq)

`--breseq_threads`

`--mapping_breseq_p` sample is not clonal. Predict polymorphic (mixed) mutations.

`--mapping_breseq_options` "-m 30 -b 20"

`--mapping_breseq_coverage` "--min-MQ 30 --min-BQ 20"

### MAPPING - MINIMAP2

see here for detailed possible options [https://lh3.github.io/minimap2/minimap2.html](https://lh3.github.io/minimap2/minimap2.html)

`--minimap2_threads`

`--mapping_minimap2_options` "--sam-hit-only --secondary=yes -ax sr"

`--mapping_minimap2_samblaster` remove duplicates with samblaster

`--mapping_minimap2_coverage` "--min-MQ 30 --min-BQ 20"

### MAPPING - BWA

see here for detailed possible options [http://bio-bwa.sourceforge.net/bwa.shtml](http://bio-bwa.sourceforge.net/bwa.shtml)

`--bwa_threads`

`--mapping_bwa_options` "-M"

`--mapping_bwa_samblaster` remove duplicates with samblaster

`--mapping_bwa_coverage` "--min-MQ 30 --min-BQ 20"

## :dango: SNPCALLING options

### SNPCALLING - 

### SNPCALLING - FREEBAYES

see here for detailed possible options [https://github.com/freebayes/freebayes](https://github.com/freebayes/freebayes)

`--freebayes_threads`

`--snpcalling_freebayes_options` "--pooled-discrete --min-alternate-fraction 0.01 --min-alternate-count 2 --min-mapping-quality 30 --min-base-quality 20"

### SNPCALLING - BCFTOOLS

see here for detailed possible options [https://samtools.github.io/bcftools/bcftools.html](https://samtools.github.io/bcftools/bcftools.html)

`--bcftools_threads`

`--snpcalling_bcftools_mpileup_options` "-C 50 -q 30 -Q 20 -d 2000 -E -a FORMAT/AD,FORMAT/DP"

`--snpcalling_bcftools_call_options` "-mAv -Ov"

`--snpcalling_bcftools_varfilter_options` "-Q 10 -d 5 -D 2000"

### SNPCALLING - GDCOMAPRE

`--gdtools_threads`

### SNPCALLING - LOFREQ

see here for detailed possible options [https://csb5.github.io/lofreq/commands/](https://csb5.github.io/lofreq/commands/)

`--lofreq_threads`

`--snpcalling_lofreq_options` "-C 5 -d 2000 -m 30 -q 20 -Q 20 -D --call-indels"

### SNPCALLING - VARSCAN

see here for detailed possible options [http://varscan.sourceforge.net/using-varscan.html](http://varscan.sourceforge.net/using-varscan.html)

`--varscan_threads`

`--snpcalling_varscan_mpileup_options` "-C 50 -q 30 -Q 20 -d 2000 -E"

`--snpcalling_varscan_snp_options` "--min-coverage 5 --min-avg-qual 20 --min-var-freq 0.01 --output-vcf 1"

`--snpcalling_varscan_indel_options` "--min-coverage 5 --min-avg-qual 20 --min-var-freq 0.01 --output-vcf 1"

## :oden: SVCALLING options

### SVCALLING - PINDEL

see here for detailed possible options [https://gmt.genome.wustl.edu/packages/pindel/user-manual.html](https://gmt.genome.wustl.edu/packages/pindel/user-manual.html)

`--pindel_threads`

`--svcalling_pindel_sam2pindel_options` "300 tag 0 Illumina-PairEnd"

`--svcalling_pindel_options` "-c ALL"

### SVCALLING - GRIDSS

see here for detailed possible options [https://github.com/PapenfussLab/gridss](https://github.com/PapenfussLab/gridss)

## ANNOTATION

### ANNOTATION - SNPEFF

see here for detailed possible options [http://pcingola.github.io/SnpEff/se_introduction/](http://pcingola.github.io/SnpEff/se_introduction/)

