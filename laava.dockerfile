# Interactive environment with scripts and extra dependencies
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

# Executable scripts
COPY bin/* /usr/local/bin
RUN chmod +x /usr/local/bin/*.py /usr/local/bin/*.R
ENV PATH "/usr/local/bin:$PATH"

WORKDIR /data/

CMD ["bash"]
