# Development environment for running the scripts, no scripts, extra dependencies
FROM --platform=linux/amd64 continuumio/miniconda3:23.10.0-1

# Set the container's timezone to match this local machine
RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone
# Silence a debconf warning
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        bedtools \
        build-essential \
        git \
        less \
        minimap2 \
        samtools \
        texlive-latex-extra \
        texlive-latex-recommended \
        vim \
    && rm -rf /var/lib/apt/lists/*

# Install directly into 'base' conda environment
COPY laava_dev.conda_env.yml ./conda_env.yml
RUN conda install -y -n base conda-libmamba-solver && conda config --set solver libmamba
RUN conda install -y -n base python=3.10
RUN conda env update -v -n base -f conda_env.yml

WORKDIR /data/

CMD ["bash"]
