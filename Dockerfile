FROM --platform=linux/amd64 continuumio/miniconda3

RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        ca-certificates \
        less \
        texlive-latex-extra \
        texlive-latex-recommended \
    && rm -rf /var/lib/apt/lists/*

# Conda environment has the name AAV.env
COPY AAV.conda_env.yml ./conda_env.yml
RUN conda install -n base conda-libmamba-solver && conda config --set solver libmamba
RUN conda env create -v -f conda_env.yml
RUN echo "conda activate AAV.env" >> ~/.bashrc

# Executable scripts for the pipeline
COPY bin/* /usr/local/bin
RUN chmod +x /usr/local/bin/*.py /usr/local/bin/*.R
ENV PATH "/usr/local/bin:$PATH"

WORKDIR /data/

CMD ["bash"]
