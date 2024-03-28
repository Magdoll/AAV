# Nextflow process environment, no scripts, just the core dependencies
FROM --platform=linux/amd64 continuumio/miniconda3:23.10.0-1

RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        samtools \
        texlive-latex-extra \
        texlive-latex-recommended \
    && rm -rf /var/lib/apt/lists/*

# Install directly into 'base' conda environment
COPY laava.conda_env.yml ./conda_env.yml
RUN conda install -y -n base conda-libmamba-solver && conda config --set solver libmamba
RUN conda install -y -n base python=3.10
RUN conda env update -v -n base -f conda_env.yml

WORKDIR /data/

CMD ["bash"]
