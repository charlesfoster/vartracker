# syntax=docker/dockerfile:1.5

FROM mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/conda/bin:$PATH

WORKDIR /app

RUN micromamba install --yes --name base \
        --channel conda-forge --channel bioconda \
        python=3.11.8 \
        snakemake=9.0.1 \
        bcftools=1.21 \
        samtools=1.21 \
        htslib=1.21 \
        fastp=0.23.4 \
        bwa=0.7.17 \
        lofreq=2.1.5 \
    && micromamba clean --all --yes

COPY --chown=mambauser:mambauser . /app

RUN pip install --no-cache-dir .

ENTRYPOINT ["vartracker"]
CMD ["--help"]
