FROM --platform=amd64 continuumio/miniconda3:latest
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y unzip

RUN conda install pandas pip -y && \
    conda install -c conda-forge google-cloud-sdk oauth2client -y && \
    pip install altocumulus firecloud tabulate

RUN mkdir /resources
COPY . /resources
ENV GOOGLE_APPLICATION_CREDENTIALS=/resources/scrnaseq-pipeline-sa.json

RUN wget https://github.com/klarman-cell-observatory/scrnaseq_pipeline/archive/master.zip -O /resources/master.zip
RUN unzip -d /resources/ /resources/master.zip
WORKDIR /resources/scrnaseq_pipeline-master/src/

