#!/bin/bash

#
# conda create -n dsub_env python=3 pip
# conda activate dsub_env
# pip install --upgrade dsub
# conda install -c conda-forge google-cloud-sdk
# gcloud auth configure-docker
#

# conda activate dsub_env

dir_name="finngen_atac"
gcp_bucket_basedir="gs://fc-secure-d4adbbf9-8265-4a5c-b14f-23a5f1b5c4f9/${dir_name}"
sample_tracking_file="${gcp_bucket_basedir}/011124_pooled_atac_count_sampletracker.csv"
project_name="finngen_atac"
email="will@broadinstitute.org"
workspace="'693-finland-v2f/Finngen'"
count_matrix_name="raw_feature_bc_matrix.h5"
steps="COUNT"
mkfastq_memory="256G"
mkfastq_diskspace="2500"
cellranger_method="broadinstitute:cumulus:Cellranger:2.2.0"
cumulus_method="broadinstitute:cumulus:cumulus:2.1.1"
cellbender_method="cellbender/remove-background/11"
cellranger_version="7.0.1"
cellranger_atac_version="2.1.0"
cellranger_arc_version="2.0.1"

current_time=$(date "+%Y.%m.%d-%H.%M.%S")

dsub --provider google-cls-v2 --project "microbiome-xavier" --regions us-east1 \
  --service-account "scrnaseq-pipeline@microbiome-xavier.iam.gserviceaccount.com" \
  --image "gcr.io/microbiome-xavier/conda-alto" --disk-size '10' --timeout '2d'\
  --logging "$gcp_bucket_basedir/logs/" \
  --command "wget https://github.com/will-broad/scrnaseq_pipeline/archive/atac.zip && unzip atac.zip && cd scrnaseq_pipeline-atac/src && python sc_pipeline.py" \
  --output PIPELINE_LOGS="$gcp_bucket_basedir/logs/execution_$current_time.log" \
  --input SAMPLE_TRACKING_FILE="$sample_tracking_file" \
  --env PROJECT_NAME="$project_name" \
  --env GCP_BUCKET_BASEDIR="$gcp_bucket_basedir" \
  --env EMAIL="$email" \
  --env TERRA_WORKSPACE="$workspace" \
  --env COUNT_MATRIX_NAME="$count_matrix_name" \
  --env STEPS="$steps" \
  --env CELLRANGER_METHOD="$cellranger_method" \
  --env CUMULUS_METHOD="$cumulus_method" \
  --env CELLBENDER_METHOD="$cellbender_method" \
  --env CELLRANGER_VERSION="$cellranger_version" \
  --env CELLRANGER_ATAC_VERSION="$cellranger_atac_version" \
  --env CELLRANGER_ARC_VERSION="$cellranger_arc_version" \
  --env MKFASTQ_DISKSPACE="$mkfastq_diskspace" \
  --env MKFASTQ_MEMORY="$mkfastq_memory"