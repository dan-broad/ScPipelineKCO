FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

# See license at https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/200017777_00_BCL_Convert_v3.10.5_Customer_Release_Notes.pdf


RUN apt-get update && apt-get install --no-install-recommends -y alien unzip build-essential dpkg-dev gnupg lsb-release procps && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk

RUN pip install --upgrade pip && \
pip install pandas

ADD bcl-convert-3.10.5-2.el7.x86_64.rpm /software/
RUN alien -i /software/bcl-convert-3.10.5-2.el7.x86_64.rpm && rm /software/bcl-convert-3*
