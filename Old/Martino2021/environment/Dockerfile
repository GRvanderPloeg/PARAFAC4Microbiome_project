# hash:sha256:29a921427c711f19a2d881bb7f2f460829a52c96968de1a34d8e11a3c42012b8
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-python3.7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

ADD qiime2-2019.7-py36-linux-conda.yml qiime2-2019.7-py36-linux-conda.yml
RUN conda env create -n qiime2-2019.7 -f qiime2-2019.7-py36-linux-conda.yml
RUN echo "conda activate qiime2-2019.7" > ~/.bashrc
ENV PATH /opt/conda/envs/qiime2-2019.7/bin:$PATH

RUN pip install -U --no-cache-dir \
    gemelli==0.0.5
RUN pip install -U --no-cache-dir \
    biopython
RUN pip install -U --no-cache-dir \
    matplotlib-venn

COPY postInstall /
RUN /postInstall
