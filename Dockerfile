# syntax=docker/dockerfile:1.5

FROM mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/conda/bin:$PATH

WORKDIR /app

RUN micromamba install --yes --name base \
        --channel conda-forge --channel bioconda \
        python=3.11 \
        snakemake \
        bcftools \
        samtools \
        tabix \
        fastp \
        bwa \
        lofreq \
    && micromamba clean --all --yes

COPY --chown=mambauser:mambauser . /app

RUN pip install --no-cache-dir .

ENTRYPOINT ["vartracker"]
CMD ["--help"]
