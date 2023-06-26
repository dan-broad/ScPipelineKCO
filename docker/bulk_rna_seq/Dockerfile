FROM gcr.io/genomics-xavier/bulk_rna_seq:latest

ADD run_DE_analysis.R /scripts/

RUN apt-get update && \
    apt-get install --no-install-recommends -y python-pip

RUN pip install pandas