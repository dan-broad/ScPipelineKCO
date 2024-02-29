FROM --platform=amd64 continuumio/miniconda3:latest
SHELL ["/bin/bash", "-c"]

# See license at https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen/200017777_00_BCL_Convert_v3.10.5_Customer_Release_Notes.pdf
ARG BCL_CONVERT_VERSION=4.2.7

RUN apt-get update && apt-get install --no-install-recommends -y alien unzip build-essential dpkg-dev gnupg lsb-release procps curl && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && apt-get install google-cloud-sdk -y
    
RUN pip install --upgrade pip && \
pip install pandas

ADD bcl-convert-${BCL_CONVERT_VERSION}-2.el7.x86_64.rpm /software/
RUN alien -i /software/bcl-convert-${BCL_CONVERT_VERSION}-2.el7.x86_64.rpm && rm /software/bcl-convert*
