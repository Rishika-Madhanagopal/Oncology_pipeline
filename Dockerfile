FROM continuumio/miniconda3

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
SHELL ["conda", "run", "-n", "scrna_env", "/bin/bash", "-c"]

WORKDIR /workspace
COPY scripts/ /workspace/scripts
