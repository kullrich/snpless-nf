## Introduction

The output of snpless-nf primarily consists of the following main components:

## Output directotry structure

The default directory structure of snpless-nf is as follows:

```
results/
├── ASSEMBLY
│   └── UNICYCLER
│       ├── <Sample1>
│       │   ├── <Sample1>.assembly.fasta
│       │   └── unicycler.log
│       └── <Sample2>
│           ├── <Sample2>.assembly.fasta
│           └── unicycler.log
├── GENMAP
│   ├── <reference>
│   ├── reference.genmap.K30.E2.bedgraph
│   ├── reference.genmap.K30.E2.chrom.sizes
│   ├── reference.genmap.K30.E2.txt
│   └── reference.genmap.K30.E2.wig
├── MAPPING
│   ├── BRESEQ
│   │   ├── <Sample1>
│   │   └── <Sample2>
│   ├── BRESEQOUT
│   ├── BWA
│   │   ├── <Sample1>
│   │   └── <Sample2>
│   ├── BWAOUT
│   │   ├── <Sample1>.bwa.mean.coverage
│   │   ├── <Sample1>.bwa.sorted.bam
│   │   ├── <Sample1>.bwa.sorted.bam.bai
│   │   ├── <Sample2>.bwa.mean.coverage
│   │   ├── <Sample2>.bwa.sorted.bam
│   │   └── <Sample2>.bwa.sorted.bam.bai
│   ├── MINIMAP2
│   │   ├── <Sample1>
│   │   └── <Sample2>
│   └── MINIMAP2OUT
│       ├── <Sample1>.minimap2.mean.coverage
│       ├── <Sample1>.minimap2.sorted.bam
│       ├── <Sample1>.minimap2.sorted.bam.bai
│       ├── <Sample2>.minimap2.mean.coverage
│       ├── <Sample2>.minimap2.sorted.bam
│       └── <Sample2>.minimap2.sorted.bam.bai
├── QC
│   ├── FASTQC
│   │   ├── <Sample1>
│   │   └── <Sample2>
│   ├── PEAR
│   │   ├── <Sample1>
│   │   └── <Sample2>
│   └── TRIM
│       ├── <Sample1>
│       └── <Sample2>
├── SNPCALLING
│   └── LOFREQ
│       ├── BRESEQ
│       │   ├── <Sample1>
│       │   ├── <Sample1>.breseq.lofreq.vcf
│       │   ├── <Sample2> 
│       │   └── <Sample2>.breseq.lofreq.vcf
│       ├── BWA
│       │   ├── <Sample1>
│       │   ├── <Sample1>.bwa.lofreq.vcf
│       │   ├── <Sample2>
│       │   └── <Sample2>.bwa.lofreq.vcf
│       └── MINIMAP2
│       │   ├── <Sample1>
│       │   ├── <Sample1>.minimap2.lofreq.vcf
│       │   ├── <Sample2>
│       │   └── <Sample2>.minimap2.lofreq.vcf
├── SVCALLING
│   └── PINDEL
│       ├── BRESEQ
│       │   ├── <Sample1>.minimap2.pindel.txt
│       │   ├── <Sample1>.minimap2.pindel_BP
│       │   ├── <Sample1>.minimap2.pindel_CloseEndMapped
│       │   ├── <Sample1>.minimap2.pindel_D
│       │   ├── <Sample1>.minimap2.pindel_INT_final
│       │   ├── <Sample1>.minimap2.pindel_INV
│       │   ├── <Sample1>.minimap2.pindel_LI
│       │   ├── <Sample1>.minimap2.pindel_RP
│       │   ├── <Sample1>.minimap2.pindel_SI
│       │   ├── <Sample1>.minimap2.pindel_TD
│       │   ├── <Sample2>.minimap2.pindel.txt
│       │   ├── <Sample2>.minimap2.pindel_BP
│       │   ├── <Sample2>.minimap2.pindel_CloseEndMapped
│       │   ├── <Sample2>.minimap2.pindel_D
│       │   ├── <Sample2>.minimap2.pindel_INT_final
│       │   ├── <Sample2>.minimap2.pindel_INV
│       │   ├── <Sample2>.minimap2.pindel_LI
│       │   ├── <Sample2>.minimap2.pindel_RP
│       │   ├── <Sample2>.minimap2.pindel_SI
│       │   └── <Sample2>.minimap2.pindel_TD
│       ├── BWA
│       │   ├── <Sample1>.minimap2.pindel.txt
│       │   ├── <Sample1>.minimap2.pindel_BP
│       │   ├── <Sample1>.minimap2.pindel_CloseEndMapped
│       │   ├── <Sample1>.minimap2.pindel_D
│       │   ├── <Sample1>.minimap2.pindel_INT_final
│       │   ├── <Sample1>.minimap2.pindel_INV
│       │   ├── <Sample1>.minimap2.pindel_LI
│       │   ├── <Sample1>.minimap2.pindel_RP
│       │   ├── <Sample1>.minimap2.pindel_SI
│       │   ├── <Sample1>.minimap2.pindel_TD
│       │   ├── <Sample2>.minimap2.pindel.txt
│       │   ├── <Sample2>.minimap2.pindel_BP
│       │   ├── <Sample2>.minimap2.pindel_CloseEndMapped
│       │   ├── <Sample2>.minimap2.pindel_D
│       │   ├── <Sample2>.minimap2.pindel_INT_final
│       │   ├── <Sample2>.minimap2.pindel_INV
│       │   ├── <Sample2>.minimap2.pindel_LI
│       │   ├── <Sample2>.minimap2.pindel_RP
│       │   ├── <Sample2>.minimap2.pindel_SI
│       │   └── <Sample2>.minimap2.pindel_TD
│       └── MINIMAP2
│           ├── <Sample1>.minimap2.pindel.txt
│           ├── <Sample1>.minimap2.pindel_BP
│           ├── <Sample1>.minimap2.pindel_CloseEndMapped
│           ├── <Sample1>.minimap2.pindel_D
│           ├── <Sample1>.minimap2.pindel_INT_final
│           ├── <Sample1>.minimap2.pindel_INV
│           ├── <Sample1>.minimap2.pindel_LI
│           ├── <Sample1>.minimap2.pindel_RP
│           ├── <Sample1>.minimap2.pindel_SI
│           ├── <Sample1>.minimap2.pindel_TD
│           ├── <Sample2>.minimap2.pindel.txt
│           ├── <Sample2>.minimap2.pindel_BP
│           ├── <Sample2>.minimap2.pindel_CloseEndMapped
│           ├── <Sample2>.minimap2.pindel_D
│           ├── <Sample2>.minimap2.pindel_INT_final
│           ├── <Sample2>.minimap2.pindel_INV
│           ├── <Sample2>.minimap2.pindel_LI
│           ├── <Sample2>.minimap2.pindel_RP
│           ├── <Sample2>.minimap2.pindel_SI
│           └── <Sample2>.minimap2.pindel_TD
├── dag.svg
├── report.html
├── timeline.html
└── trace.txt
work/
```

The default name of the output directory (unless otherwise specified) will be `./results/` as specified with the `--output` option.

## Primary Output Directories

## Secondary Output Directories


