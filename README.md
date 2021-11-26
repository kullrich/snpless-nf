# snpless-nf
snpless-nf - A Nextflow pipeline for time-course analysis with bacterial NGS whole-genome data.

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-brightgreen.svg)](http://nextflow.io)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![MITlicense](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

## Introduction

## Pipeline summary

1. QC
    1. FASTQC [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    2. TRIM [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
    3. PEAR [pear](https://cme.h-its.org/exelixis/web/software/pear/)
2. GENMAP [GenMap](https://github.com/cpockrandt/genmap)
3. ASSEMBLY
    1. UNICYCLER [Unicycler](https://github.com/rrwick/Unicycler)
    2. PROKKA [prokka]((https://github.com/tseemann/prokka)
4. MAPPING
    1. BRESEQ [breseq](https://github.com/barricklab/breseq) >> SAMTOOLS [samtools](https://github.com/samtools/samtools) add read group
    2. MINIMAP2 [minimap2](https://github.com/lh3/minimap2) >> SAMBLASTER [samblaster](https://github.com/GregoryFaust/samblaster) remove duplicates
    3. BWA [BWA](http://bio-bwa.sourceforge.net/) >> SAMBLASTER [samblaster](https://github.com/GregoryFaust/samblaster) remove duplicates
    4. COVERAGE [samtools](https://github.com/samtools/samtools)
5. SNPCALLING
    1. MPILEUP [samtools](https://github.com/samtools/samtools)
    2. BCFTOOLS [bcftools](https://github.com/samtools/bcftools)
    3. FREEBAYES [freebayes](https://github.com/freebayes/freebayes)
    4. LOFREQ [LoFreq](http://csb5.github.io/lofreq/)
    5. VARSCAN [varscan](http://dkoboldt.github.io/varscan/)
6. SVCALLING
    1. PINDEL [pindel](https://github.com/genome/pindel)
    2. BREAKDANCER [BreakDancer](https://github.com/genome/breakdancer)
7. FILTERING/MERGING
    1. GDCOMPARE [gdtools](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_usage.html)
    2. BEDTOOLS [bedtools](https://bedtools.readthedocs.io/en/latest/)
8. ANNOTATION
    1. SNPEFF [SnpEff](http://pcingola.github.io/SnpEff/)
9. PLOTTING
    1. PLOT [R](https://cran.r-project.org/)

Addtional Tools used for data conversion:

- HTSLIB [htslib]((https://github.com/samtools/htslib)
- 

## Quickstart

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) (>=21.10.0)

Install Nextflow by using the following command:

```
curl -s https://get.nextflow.io | bash
```

or

Install Nextflow by using conda:

```
conda create -n nf python=3
conda activate nf
conda install -c bioconda nextflow
```

2. Download the pipeline

```
git clone https://github.com/kullrich/snpless-nf.git
```

3. Test the pipeline on an minimal dataset with a single command:

Using <nextflow> conda environment:

```
conda activate nf
nextflow run snpless-nf -profile test
```

4. Start running your own analysis:

Check the necessary input files!

```
nextflow run snpless-nf --input <samples.tsv> --reference <genome.fna> --gff3 <genome.gff3> --proteins <genome.gbff>
```

## Full example dataset

### Get example files (8.6 GB)

Download via wget:

```
cd snpless-nf/examples
wget 
tar -xvf behringer2018.tar.gz
```

Download via weblink:

[behringer2018 - samples 113, 129, 221](https://ftp.evolbio.mpg.de/main.html?download&weblink=74b3a1f98426435d16a97bcc8e55b400&realfilename=behringer2018.tar.gz)

### Run full example dataset

Using <nextflow> conda environment:

```
conda activate nf
nextflow run snpless-nf --input behringer2018/behringer2018_113.txt --reference behringer2018/GCF_000005845.2_ASM584v2_genomic.fna --gff3 GCF_000005845.2_ASM584v2_genomic.gff --proteins behringer2018/GCF_000005845.2_ASM584v2_genomic.gbff
```

## Pipeline usage

see a detailed description here:

### Input files

## Pipeline parameters

see a detailed description here:

## Pipeline output

see a detailed description here:

## Licence

MIT (see LICENSE)

## Contributing Code

If you would like to contribute to snpless-nf, please file an issue so that one can establish a statement of need, avoid redundant work, and track progress on your contribution.

Before you do a pull request, you should always file an issue and make sure that someone from the snpless-nf developer team agrees that it’s a problem, and is happy with your basic proposal for fixing it.

Once an issue has been filed and we've identified how to best orient your contribution with package development as a whole, [fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [main repo](https://github.com/kullrich/snpless-nf.git), branch off a [feature branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches) from `master`, [commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/committing-and-reviewing-changes-to-your-project) and [push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository) your changes to your fork and submit a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests) for `snpless-nf:master`.

By contributing to this project, you agree to abide by the Code of Conduct terms.

## Bug reports

Please report any errors or requests regarding [snpless-nf](https://github.com/kullrich/snpless-nf) to Kristian Ullrich (ullrich@evolbio.mpg.de)

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://github.com/kullrich/snpless-nf/-/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

## References

Behringer, Megan G., et al. "Escherichia coli cultures maintain stable subpopulation structure during long-term evolution." Proceedings of the National Academy of Sciences 115.20 (2018): E4642-E4650. [https://www.pnas.org/content/115/20/E4642.short](https://www.pnas.org/content/115/20/E4642.short)

