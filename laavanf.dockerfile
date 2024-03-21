# Nextflow process environment, no scripts, just the core dependencies
FROM --platform=linux/amd64 continuumio/miniconda3

RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        ca-certificates \
        less \
        texlive-latex-extra \
        texlive-latex-recommended \
    && rm -rf /var/lib/apt/lists/*

# Conda environment has the name 'laava'
COPY laava.conda_env.yml ./conda_env.yml
RUN conda install -y -n base conda-libmamba-solver && conda config --set solver libmamba
RUN conda env create -y -v -f conda_env.yml
RUN echo "conda activate laava" >> ~/.bashrc

WORKDIR /data/

CMD ["bash"]
