FROM nfcore/base
LABEL authors="Kristian K Ullrich" \
      email="ullrich@evolbio.mpg.de" \
      description="Docker image containing all requirements for the snpless-nf pipeline"

COPY qc.yml /
     mappability.yml /
     assembly_unicycler.yml /
     assembly_prokka.yml /
     mapping_minimap2.yml /
     mapping_bwa.yml /
     mapping_breseq.yml /
     mapping_samblaster.yml /
     snpcalling_bcftools.yml /
     snpcalling_varscan.yml /
     snpcalling_lofreq.yml /
     svcalling_pindel.yml /
     svcalling_breakdancer.yml /
     svcalling_delly.yml /
     svcalling_cnvator.yml /
RUN conda env create -f /qc.yml && conda clean -a && \
    conda env create -f /mappability.yml && conda clean -a && \
    conda env create -f /assembly_unicycler.yml && conda clean -a && \
    conda env create -f /assembly_prokka.yml && conda clean -a && \
    conda env create -f /mapping_minimap2.yml && conda clean -a && \
    conda env create -f /mapping_bwa.yml && conda clean -a && \
    conda env create -f /mapping_samblaster.yml && conda clean -a && \
    conda env create -f /snpcalling_bcftools.yml && conda clean -a && \
    conda env create -f /snpcalling_varscan.yml && conda clean -a && \
    conda env create -f /snpcalling_lofreq.yml && conda clean -a && \
    conda env create -f /svcalling_pindel.yml && conda clean -a && \
    conda env create -f /svcalling_breakdancer.yml && conda clean -a && \
    conda env create -f /svcalling_delly.yml && conda clean -a && \
    conda env create -f /svcalling_cnvator.yml && conda clean -a
ENV PATH $PATH
