import datetime
import logging
import concurrent.futures
import threading
import os
import pandas as pd
import steps
from consts import MULTIOME, RNA, ATAC
from utils import *

"""
Config Section - Modify this section only
"""
project_name = os.getenv("PROJECT_NAME", default="Gut_eQTL")
sample_tracking_file = os.getenv("SAMPLE_TRACKING_FILE", default="../data/sampletracking_multiome.csv")
gcp_basedir = os.getenv("GCP_BUCKET_BASEDIR", default="gs://fc-secure-1620151c-e00c-456d-9daf-4d222e1cab18/Gut_eQTL")
email = os.getenv("EMAIL", default="dchafamo@broadinstitute.org")
alto_workspace = os.getenv("TERRA_WORKSPACE", default="'kco-tech/Gut_eQTL'")
count_matrix_name = os.getenv("COUNT_MATRIX_NAME", default="filtered_feature_bc_matrix.h5")
steps_to_run = os.getenv("STEPS", default="BCL_CONVERT,COUNT,CUMULUS").split(',')
mkfastq_disk_space = int(os.getenv("MKFASTQ_DISKSPACE", default=1500))
mkfastq_memory = os.getenv("MKFASTQ_MEMORY", default="120G")
cellbender_method = os.getenv("CELLBENDER_METHOD", default="cellbender/remove-background/13")
cellbender_version = os.getenv("CELLBENDER_VERSION", default="0.3.0")
cumulus_method = os.getenv("CUMULUS_METHOD", default="broadinstitute:cumulus:cumulus:2.1.1")
cellranger_method = os.getenv("CELLRANGER_METHOD", default="broadinstitute:cumulus:Cellranger:2.1.1")
cellranger_version = os.getenv("CELLRANGER_VERSION", default="7.0.1")
cellranger_atac_version = os.getenv("CELLRANGER_ATAC_VERSION", default="2.1.0")
cellranger_arc_version = os.getenv("CELLRANGER_ARC_VERSION", default="2.0.1")
# BCL Convert configs
bcl_convert_method = os.getenv("BCL_CONVERT_METHOD", default="kco/bcl_convert/11")
bcl_convert_version = os.getenv("BCL_CONVERT_VERSION", default="4.2.7")
bcl_convert_disk_space = int(os.getenv("BCL_CONVERT_DISK_SPACE", default="1500"))
bcl_convert_memory = int(os.getenv("BCL_CONVERT_MEMORY", default="120"))
bcl_convert_cpu = int(os.getenv("BCL_CONVERT_NUM_CPU", default="32"))
bcl_convert_strict_mode = eval(os.getenv("BCL_CONVERT_STRICT_MODE", default="False"))
bcl_convert_file_format_version = os.getenv("BCL_CONVERT_FILE_FORMAT_VERSION", default="2")
bcl_convert_lane_splitting = eval(os.getenv("BCL_CONVERT_LANE_SPLITTING", default="False"))
bcl_convert_num_lanes = int(os.getenv("NUM_LANES_FLOWCELL", default="0"))
bcl_convert_gex_i5_index_key = os.getenv("GEX_I5_INDEX_KEY", default='index2_workflow_a(i5)')
bcl_convert_docker_registry = os.getenv("BCL_CONVERT_DOCKER_REGISTRY", default="us-docker.pkg.dev/microbiome-xavier/broad-microbiome-xavier")
# global configs
terra_timeout = int(os.getenv("TERRA_TIMEOUT", default='18000'))

"""
Set global variables
"""
max_parallel_threads = 50
cellbender_matrix_name = "out_FPR_0.01_filtered.h5"
cwd = os.getcwd()
basedir = cwd + "/" + project_name + "/sc_processed"
os.makedirs(basedir, exist_ok=True)
directories = build_directories(basedir)

"""
Preprocess Sample tracking file and Sanity check columns
"""

master_tracking = pd.read_csv(sample_tracking_file)
master_tracking['seq_dir'] = master_tracking['seq_dir'].apply(lambda sd: sd[:-1] if sd.endswith('/') else sd)
master_tracking['Sample'] = master_tracking['sampleid']
project = master_tracking[master_tracking.run_pipeline]['project'].tolist()[0]
buckets = build_buckets(gcp_basedir, project)
alto_dirs = build_alto_folders(buckets)
log_file = os.getenv("PIPELINE_LOGS", default='{}/{}.log'.format(basedir, project_name))

sample_sheet_columns = [
    'date', 'run_pipeline', 'Channel Name', 'Sample', 'sampleid', 'method', 'sub_method',
    'condition', 'replicate', 'tissue', 'Lane', 'Index', 'instrument_platform', 'instrument_type',
    'create_fastq_for_index_reads', 'trim_umi', 'override_cycles', 'project', 'reference',
    'introns', 'chemistry', 'flowcell', 'seq_dir', 'min_umis', 'min_genes', 'percent_mito', 
    'cellbender_expected_cells', 'cellbender_total_droplets_included', 'cellbender_learning_rate',
    'cellbender_force_cell_umi_prior', 'cellbender_force_empty_umi_prior'
]

for col in sample_sheet_columns:
    if col not in master_tracking.columns:
        logging.error(f"Missing columns: {col} in samplesheet. Exiting.")
        exit(1)


def process_bcl_convert(sample_tracking):
    sample_tracking = sample_tracking.reset_index(drop=True) #reset index for multiome 
    threading.current_thread().name = f'Thread: bcl_convert ({sample_tracking["sub_method"][0]})'

    env_vars = {
           "software_version": bcl_convert_version, 
           "delete_input_dir": False, 
           "disk_space": bcl_convert_disk_space, 
           "memory": bcl_convert_memory, 
           "cpu": bcl_convert_cpu, 
           "strict_mode": bcl_convert_strict_mode, 
           "file_format_version": bcl_convert_file_format_version,
           "no_lane_splitting": not bcl_convert_lane_splitting,
           "num_lanes": bcl_convert_num_lanes,
           "gex_i5_index_key": bcl_convert_gex_i5_index_key,
           "docker_registry": bcl_convert_docker_registry
       }

    paths = steps.upload_bcl_convert_input(
        sample_tracking, 
        buckets, 
        directories,
        env_vars
    )
    
    steps.run_bcl_convert(
        directories, 
        buckets, 
        paths, 
        bcl_convert_method, 
        alto_workspace,
        terra_timeout
    )
    
    if env_vars["no_lane_splitting"]:
        steps.move_fastqs_to_sample_dir(directories, buckets, sample_tracking)

def process_rna_flowcell(seq_dir):
    """
    Initiate pipeline for a set of samples within a single Flowcell.
    :param seq_dir: GCP Cloud Storage link to raw BCL directory
    """
    sample_tracking = master_tracking[master_tracking.run_pipeline &
                                      (master_tracking.seq_dir == seq_dir)]
    
    threading.current_thread().name = 'Thread:' + sample_tracking['flowcell'].iloc[0]
    logging.info("Started processing samples in {}".format(seq_dir))

    sample_tracking = sample_tracking[sample_sheet_columns]

    sample_dicts = build_sample_dicts(sample_tracking, sample_tracking['sampleid'].tolist())

    if "BCL_CONVERT" in steps_to_run:
        process_bcl_convert(sample_tracking)

    # DEPRECATED
    if "MKFASTQ" in steps_to_run:

        steps.upload_cellranger_mkfastq_input(
            buckets,
            directories,
            sample_tracking,
            cellranger_version,
            cellranger_atac_version,
            mkfastq_disk_space,
            mkfastq_memory
        )

        steps.run_cellranger_mkfastq(
            directories,
            sample_tracking,
            alto_workspace,
            cellranger_method,
            alto_dirs['alto_fastqs'],
            terra_timeout
        )

    if "COUNT" in steps_to_run:

        steps.upload_cellranger_count_input(
            buckets,
            directories,
            sample_dicts,
            sample_tracking,
            cellranger_version,
            cellranger_atac_version
        )

        steps.run_cellranger_count(
            directories,
            sample_dicts,
            sample_tracking,
            alto_workspace,
            cellranger_method,
            alto_dirs['alto_counts'],
            terra_timeout
        )

    if "CUMULUS" in steps_to_run:

        steps.upload_cumulus_samplesheet(
            buckets,
            directories,
            sample_dicts,
            sample_tracking,
            count_matrix_name
        )

        steps.run_cumulus(
            directories,
            sample_dicts,
            sample_tracking,
            alto_workspace,
            cumulus_method,
            alto_dirs['alto_results'],
            terra_timeout
        )

    if "CELLBENDER" in steps_to_run:

        steps.upload_cell_bender_input(
            buckets,
            directories,
            sample_dicts,
            sample_tracking,
            count_matrix_name,
            cellbender_version
        )

        steps.run_cellbender(
            directories,
            sample_dicts,
            sample_tracking,
            alto_workspace,
            cellbender_method,
            alto_dirs['alto_cellbender'],
            terra_timeout
        )

    if "CELLBENDER_CUMULUS" in steps_to_run:

        steps.upload_post_cellbender_cumulus_input(
            buckets,
            directories,
            sample_dicts,
            sample_tracking,
            cellbender_matrix_name
        )

        steps.run_cumulus_post_cellbender(
            directories,
            sample_dicts,
            sample_tracking,
            alto_workspace,
            cumulus_method,
            alto_dirs['alto_results'],
            terra_timeout
        )


def process_multiome():
    """
    Initiate pipeline for all multiome assay samples.
    """
    sample_tracking = master_tracking[master_tracking.run_pipeline &
                                      (master_tracking.method == MULTIOME)]

    threading.current_thread().name = 'Thread: MULTIOME'
    logging.info("Started processing multiome samples")

    sample_tracking = sample_tracking[sample_sheet_columns]

    if "BCL_CONVERT" in steps_to_run:
        # FASTQ generation can be parallelized 
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel_threads) as executor:
            sub_methods = set(sample_tracking['sub_method'])
            futures = []
            for m in sub_methods:
                futures.append(executor.submit(process_bcl_convert, sample_tracking[sample_tracking.sub_method == m]))

            for future in concurrent.futures.as_completed(futures):
                try:
                    logging.info(future.result())
                except Exception as e:
                    logging.error(e)

    # add path to fastq for each sample before running cellranger 
    get_run_id = lambda sample: os.path.basename(sample['seq_dir'])
    get_fastq_path = lambda r: f"{buckets['fastqs']}/{r['sub_method']}/{get_run_id(r)}_fastqs/sample_fastqs/{r['sampleid']}"
    sample_tracking['fastq_dir'] = sample_tracking.apply(get_fastq_path, axis=1)

    steps.upload_cellranger_arc_samplesheet(buckets, directories, sample_tracking, cellranger_arc_version,
                                            mkfastq_disk_space, mkfastq_memory, steps_to_run)
    steps.run_cellranger_arc(buckets, directories, cellranger_method, alto_workspace, terra_timeout)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)-7s | %(threadName)-15s | %(levelname)-5s | %(message)s",
                        level=logging.INFO, datefmt="%m-%d %H:%M", filename=log_file, filemode='w')

    start_time = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    logging.info("Running scRNASeq pipeline for project {} on {}".format(project, start_time))
    logging.info("GCP User: {}".format(email))
    logging.info("GCP bucket dir: {}".format(gcp_basedir))
    logging.info("Workspace: {}".format(alto_workspace))
    logging.info("Count matrix name: {}".format(count_matrix_name))
    logging.info("Steps: {}".format(steps_to_run))
    logging.info("Master sample tracking file: \n\n {} \n".format(master_tracking.to_markdown()))

    method = set(master_tracking[master_tracking.run_pipeline]['method'])
    logging.info(f'Methods = {method}')
    if RNA in method or ATAC in method:
        logging.info('Processing RNA Seq and ATAC Seq Samples.')
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel_threads) as executor:
            seq_dirs = set(master_tracking[master_tracking.run_pipeline & ((master_tracking.method == RNA) | (master_tracking.method == ATAC))]['seq_dir'])
            results = executor.map(process_rna_flowcell, seq_dirs)
            for res in results:
                try:
                    logging.info(res.result())
                except Exception as e: 
                    logging.error(e)
    if MULTIOME in method:
        logging.info('Processing Multiome Samples.')
        process_multiome()
