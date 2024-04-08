import datetime
import logging

from utils import *
import pandas as pd
import concurrent.futures
import steps
import threading
import os

"""
Config Section - Modify this section only
"""
project_name = os.getenv("PROJECT_NAME", default="Gut_eQTL")
sample_tracking_file = os.getenv("SAMPLE_TRACKING_FILE", default="../data/sampletracking_multiome.csv")
gcp_basedir = os.getenv("GCP_BUCKET_BASEDIR", default="gs://fc-secure-1620151c-e00c-456d-9daf-4d222e1cab18/Gut_eQTL")
email = os.getenv("EMAIL", default="will@broadinstitute.org")
alto_workspace = os.getenv("TERRA_WORKSPACE", default="'kco-tech/Gut_eQTL'")
count_matrix_name = os.getenv("COUNT_MATRIX_NAME", default="filtered_feature_bc_matrix.h5")
steps_to_run = os.getenv("STEPS", default="CELLBENDER").split(',')
mkfastq_disk_space = int(os.getenv("MKFASTQ_DISKSPACE", default=1500))
mkfastq_memory = os.getenv("MKFASTQ_MEMORY", default="120G")
cellbender_method = os.getenv("CELLBENDER_METHOD", default="broadinstitute:cumulus:CellBender:2.5.0")
cumulus_method = os.getenv("CUMULUS_METHOD", default="broadinstitute:cumulus:cumulus:2.1.1")
cellranger_method = os.getenv("CELLRANGER_METHOD", default="broadinstitute:cumulus:Cellranger:2.1.1")
cellranger_version = os.getenv("CELLRANGER_VERSION", default="7.0.1")
cellranger_atac_version = os.getenv("CELLRANGER_ATAC_VERSION", default="2.1.0")
cellranger_arc_version = os.getenv("CELLRANGER_ARC_VERSION", default="2.0.1")
"""
Set global variables
"""
max_parallel_threads = 50
cellbender_matrix_name = "out_FPR_0.01_filtered.h5"
cwd = os.getcwd()
basedir = cwd + "/" + project_name + "/sc_processed"
os.makedirs(basedir, exist_ok=True)
directories = build_directories(basedir)
MULTIOME = 'multiome'
RNA = 'rna'
ATAC = 'atac'

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
    'date', 'run_pipeline', 'Channel Name', 'Sample', 'sampleid', 'method', 'sub_method', 'condition',
    'replicate', 'tissue', 'Lane', 'Index', 'project', 'reference', 'introns', 'chemistry', 'flowcell',
    'seq_dir', 'min_umis', 'min_genes', 'percent_mito', 'cellbender_expected_cells',
    'cellbender_total_droplets_included'
]

for col in sample_sheet_columns:
    if col not in master_tracking.columns:
        logging.error(f"Missing columns: {col} in samplesheet. Exiting.")
        exit(1)


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
            alto_dirs['alto_fastqs']
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
            alto_dirs['alto_counts']
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
            alto_dirs['alto_results']
        )

    if "CELLBENDER" in steps_to_run:

        steps.upload_cell_bender_input(
            buckets,
            directories,
            sample_dicts,
            sample_tracking,
            count_matrix_name
        )

        steps.run_cellbender(
            directories,
            sample_dicts,
            sample_tracking,
            alto_workspace,
            cellbender_method,
            alto_dirs['alto_cellbender']
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
            alto_dirs['alto_results']
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

    steps.upload_cellranger_arc_samplesheet(buckets, directories, sample_tracking, cellranger_arc_version,
                                            mkfastq_disk_space, mkfastq_memory, steps_to_run)
    steps.run_cellranger_arc(buckets, directories, cellranger_method, alto_workspace)


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
            executor.map(process_rna_flowcell, seq_dirs)
    if MULTIOME in method:
        logging.info('Processing Multiome Samples.')
        process_multiome()

