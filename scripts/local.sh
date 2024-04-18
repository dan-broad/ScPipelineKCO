#!/bin/bash

#
# conda create -n dsub_env python=3 pip
# conda activate dsub_env
# pip install --upgrade dsub
# conda install -c conda-forge google-cloud-sdk
# gcloud auth configure-docker
#

# conda activate dsub_env

dir_name="tutorial-rna"
gcp_bucket_basedir="gs://fc-secure-15bf93cd-d43c-4a70-b7de-0ee36bf3a52a/${dir_name}"
sample_tracking_file="${gcp_bucket_basedir}/test_pipeline.csv"
project_name="test_pipeline"
email="dsrirang@broadinstitute.org"
workspace="'kco-tech/sc_pipeline_tutorial'"
count_matrix_name="raw_feature_bc_matrix.h5"
steps="BCL_CONVERT,COUNT,CELLBENDER,CELLBENDER_CUMULUS"
mkfastq_memory="120G"
mkfastq_diskspace="1500"
cellranger_method="broadinstitute:cumulus:Cellranger:2.1.1"
cumulus_method="broadinstitute:cumulus:cumulus:2.1.1"
cellbender_method="cellbender/remove-background/13"
cellranger_version="7.0.1"
cellranger_atac_version="2.1.0"
cellranger_arc_version="2.0.1"
bcl_convert_method="kco/bcl_convert/12"
bcl_convert_version="4.2.7"
bcl_convert_disk_space="1500"
bcl_convert_cpu="32"
bcl_convert_strict_mode=False
bcl_convert_file_format_version="2"
bcl_convert_memory="120"
bcl_convert_lane_splitting=False
bcl_convert_docker_registry="us-docker.pkg.dev/microbiome-xavier/broad-microbiome-xavier"
bcl_convert_num_lanes_flowcell="0" # Optional: only needed when using * in sample sheet for lanes and no_lane_splitting == false
gex_i5_index_key="index2_workflow_b(i5)"

current_time=$(date "+%Y.%m.%d-%H.%M.%S")

dsub --provider google-cls-v2 --project "microbiome-xavier" --regions us-east1 \
  --service-account "scrnaseq-pipeline@microbiome-xavier.iam.gserviceaccount.com" \
  --image "gcr.io/microbiome-xavier/conda-alto" --disk-size '10' --boot-disk-size '30' --timeout '2d'\
  --logging "$gcp_bucket_basedir/logs/" \
  --command "wget https://github.com/klarman-cell-observatory/scrnaseq_pipeline/archive/master.zip && unzip master.zip && cd scrnaseq_pipeline-master && python src/sc_pipeline.py" \
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
  --env MKFASTQ_MEMORY="$mkfastq_memory" \
  --env BCL_CONVERT_METHOD="$bcl_convert_method" \
  --env BCL_CONVERT_VERSION="$bcl_convert_version" \
  --env BCL_CONVERT_DISK_SPACE="$bcl_convert_disk_space" \
  --env BCL_CONVERT_MEMORY="$bcl_convert_memory" \
  --env BCL_CONVERT_NUM_CPU="$bcl_convert_cpu" \
  --env BCL_CONVERT_STRICT_MODE="$bcl_convert_strict_mode" \
  --env BCL_CONVERT_FILE_FORMAT_VERSION="$bcl_convert_file_format_version" \
  --env NUM_LANES_FLOWCELL="$bcl_convert_num_lanes_flowcell" \
  --env BCL_CONVERT_DOCKER_REGISTRY="$bcl_convert_docker_registry" \
  --env BCL_CONVERT_LANE_SPLITTING="$bcl_convert_lane_splitting" \
  --env GEX_I5_INDEX_KEY="$gex_i5_index_key"