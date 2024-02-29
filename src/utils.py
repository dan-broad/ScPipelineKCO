import logging
import os
import re
import time
import sys
import firecloud.api as fapi
import subprocess
import threading
import json
import pandas as pd 
from string import Template
import re
from consts import *

alto_lock = threading.Lock()


def build_directories(basedir):
    directories = {
        'scripts': basedir + "/scripts",
        'fastqs': basedir + "/fastqs",
        'counts': basedir + "/counts",
        'results': basedir + "/cumulus",
        'cellranger_arc': basedir + "/cellranger_arc",
        'cellbender': basedir + "/cellbenderV2",
        'cellbender_results': basedir + "/cellbenderV2_cumulus",
        'bcl_convert': basedir + "/bcl_convert" 
    }
    for directory in directories.values():
        if not os.path.exists(directory):
            os.makedirs(directory)
    return directories


def build_buckets(gcp_basedir, project):
    return {
        'fastqs': gcp_basedir + "/fastqs_" + project,
        'counts': gcp_basedir + "/counts_" + project,
        'results': gcp_basedir + "/cumulus_" + project,
        'cellranger_arc': gcp_basedir + "/cellranger_arc_" + project,
        'cellbender': gcp_basedir + "/cellbenderv2_" + project,
        'cellbender_results': gcp_basedir + "/cellbenderv2_cumulus_" + project,
        'bcl_convert': gcp_basedir + "/bcl_convert_" + project
    }


def build_alto_folders(buckets):
    return {
        'alto_fastqs': re.sub(r'^gs://.*/', "", buckets['fastqs']),
        'alto_counts': re.sub(r'^gs://.*/', "", buckets['counts']),
        'alto_results': re.sub(r'^gs://.*/', "", buckets['results']),
        'alto_cellranger_arc': re.sub(r'^gs://.*/', "", buckets['cellranger_arc']),
        'alto_cellbender': re.sub(r'^gs://.*/', "", buckets['cellbender']),
        'alto_cellbender_results': re.sub(r'^gs://.*/', "", buckets['cellbender_results'])
    }


def build_sample_dicts(sample_tracking, sampleids):
    sample_dict = dict([(sample, []) for sample in sampleids])
    mkfastq_dict = dict()
    cumulus_dict = dict()
    cellbender_dict = dict()
    cellranger_dict = dict()
    for index, row in sample_tracking.iterrows():
        sample_dict[row['sampleid']].append(row['Sample'])
        mkfastq_dict[row['Sample']] = [row['Lane'], row['Index'], row['reference'], row['chemistry'], row['method']]
        cumulus_dict[row['sampleid']] = [row['min_umis'], row['min_genes'], row['percent_mito']]
        cellbender_dict[row['sampleid']] = [row['cellbender_expected_cells'], row['cellbender_total_droplets_included']]
        cellranger_dict[row['sampleid']] = [row['introns']]

    return {
        'sample': sample_dict,
        'mkfastq': mkfastq_dict,
        'cumulus': cumulus_dict,
        'cellbender': cellbender_dict,
        'cellranger': cellranger_dict
    }


def execute_alto_command(run_alto_file):
    with alto_lock:
        command = "bash %s" % run_alto_file
        logging.info("Executing command: `{}`".format(command))
        with open(run_alto_file, 'r') as f:
            logging.info(f"Alto file:\n {f.read()}")
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
        alto_outputs = [status_url for status_url in result.stdout.decode('utf-8').split("\n") if "http" in status_url]

    if len(alto_outputs) == 0:
        logging.info("Alto submission status url not found. %s" % result)
        sys.exit()

    for status_url in alto_outputs:
        wait_for_terra_submission(status_url)


def wait_for_terra_submission(status_url):
    logging.info("Job status: %s" % status_url)
    entries = status_url.split('/')
    workspace_namespace, workspace_name, submission_id = [entries[idx] for idx in [-4, -3, -1]]
    response = fapi.get_submission(workspace_namespace, workspace_name, submission_id)
    log_workflow_details(response)
    start_time = time.time()
    while response.json()['status'] != 'Done':
        status = [v for k, v in response.json().items() if k in ['status', 'submissionId']]
        logging.info("Job status: %s " % status)
        time.sleep(TERRA_POLL_SPACER)
        response = fapi.get_submission(workspace_namespace, workspace_name, submission_id)
        if (time.time() - start_time) > TERRA_TIMEOUT:
            logging.info("Terra pipeline took too long to complete.")
            sys.exit()
    status = {k: v for k, v in response.json().items() if k in ['status', 'submissionDate', 'submissionId']}
    logging.info("Job status: %s \n" % status)

    for workflow in response.json()['workflows']:
        if workflow['status'] != 'Succeeded':
            logging.info("Terra pipeline failed.")
            sys.exit()
    logging.info("Terra job complete: %s" % status)


def bash_execute_file(file):
    command = "bash %s" % file
    logging.info("Executing command: `{}`".format(command))
    with open(file, 'r') as f:
        logging.info(f"Bash file:\n {f.read()}")
    subprocess.run(command, shell=True, stdout=sys.stdout, stderr=sys.stderr, check=True)


def log_workflow_details(response):
    try:
        formatted_response = response.json()
        formatted_response['workflow_details'] = []
        for workflow in response.json()['workflows']:
            workflow['input'] = {}
            for res in workflow['inputResolutions']:
                workflow['input'][res['inputName']] = res['value']
            del workflow['inputResolutions']
            formatted_response['workflow_details'].append(workflow)
        del formatted_response['workflows']
        for line in json.dumps(formatted_response, sort_keys=True, indent=4).split('\n'):
            logging.info(line)
    except:
        logging.info('Unable to gather workflow input.')

def create_bcl_convert_sample_sheet(path, sub_method, env_vars, sample_tracking):
    if sub_method == 'atac':
        sample_sheet = get_atac_sample_sheet(env_vars, sample_tracking, env_vars['num_lanes'])
    elif sub_method == 'rna':
        sample_sheet = get_rna_sample_sheet(env_vars, sample_tracking, env_vars['num_lanes'])
    
    with open(path, 'w') as f:
        f.write(sample_sheet)

def get_rna_sample_sheet(env_vars, sample_tracking, num_lanes):
    dir = os.path.dirname(os.path.realpath(__file__))
    with open(f'{dir}/templates/bcl_convert_gex_sample_sheet_template.csv') as f:
        template = Template(f.read())
        sample_sheet = template.safe_substitute(env_vars)
    
    columns = ["sampleid", "Index", "Index2"]
    columns = columns + ["Lane"] if not env_vars.get("no_lane_splitting") else columns
    samples_with_indices = sample_tracking[columns]

    if not env_vars.get("no_lane_splitting"):
        samples_with_indices = apply_lane_splits(samples_with_indices, num_lanes)

    samples_with_indices = samples_with_indices.rename(columns={"sampleid": "Sample_ID", "Index": "index", "Index2": "index2"})
    return sample_sheet + samples_with_indices.to_csv(index=False,)

def get_atac_sample_sheet(env_vars, sample_tracking, num_lanes):
    dir = os.path.dirname(os.path.realpath(__file__))
    with open(f'{dir}/templates/bcl_convert_atac_sample_sheet_template.csv') as f:
        template = Template(f.read())
        sample_sheet = template.safe_substitute(env_vars)

    columns = ["Lane", "sampleid", "Index", "Index2", "Index3", "Index4"]
    samples_with_indices = sample_tracking.get(columns)
    flattend_sample_indices = flatten_sample_indices(samples_with_indices)

    if env_vars.get("no_lane_splitting"): 
        flattend_sample_indices.drop(columns='Lane', inplace=True)
    else:
        flattend_sample_indices = apply_lane_splits(flattend_sample_indices, num_lanes)

    return sample_sheet + flattend_sample_indices.to_csv(index=False)

def flatten_sample_indices(sample_tracking):
    flattened = pd.DataFrame(columns=["Lane", "Sample_ID", "index"])
    for _, r in sample_tracking.iterrows():
        sampleID = r['sampleid']
        sample_df = pd.DataFrame([
            [r['Lane'], sampleID, r['Index']],
            [r['Lane'], sampleID, r['Index2']],
            [r['Lane'], sampleID, r['Index3']],
            [r['Lane'], sampleID, r['Index4']],
        ], columns=["Lane", "Sample_ID", "index"])
        flattened = pd.concat([flattened, sample_df], ignore_index=True)
    return flattened

def apply_lane_splits(sample_tracking, num_lanes):
    lane_split = []
    
    for _, sample in sample_tracking.iterrows():
        laneVal = str(sample['Lane'])
        drop_lane = sample.drop(labels=['Lane'])
        if len(laneVal) == 1 and laneVal.isdecimal(): # no need to split
            lane_split.append([sample['Lane']] + drop_lane.to_list())
        elif '-' in laneVal:
            start, end = laneVal.split('-')
            r = range(int(start), int(end)+1)
            for i in r:
                lane_split.append([i] + drop_lane.to_list())
        elif laneVal == '*':
            for i in range(num_lanes):
                lane_split.append([i+1] + drop_lane.to_list())
        
    return pd.DataFrame(lane_split, columns=['Lane'] + sample_tracking.columns.drop('Lane').to_list())

def create_bcl_convert_params(file_path, env_vars, input_dir, fastq_output_dir, sample_sheet_path):
    with open(file_path, "w") as f: 
        params = {
                "bclconvert.bcl_convert_version": f"{env_vars['software_version']}",
                "bclconvert.delete_input_bcl_directory": env_vars['delete_input_dir'],
                "bclconvert.disk_space": env_vars['disk_space'],
                "bclconvert.docker_registry": env_vars['docker_registry'],
                "bclconvert.input_bcl_directory": input_dir,
                "bclconvert.memory": f"{env_vars['memory']}G",
                "bclconvert.no_lane_splitting": env_vars["no_lane_splitting"],
                "bclconvert.num_cpu": env_vars['cpu'],
                "bclconvert.output_directory": fastq_output_dir,
                # "bclconvert.run_bcl_convert.run_id": "${}",
                "bclconvert.sample_sheet": sample_sheet_path,
                "bclconvert.strict_mode": env_vars['strict_mode']
                # "bclconvert.zones": "${}"
            }
        
        contents = json.dumps(params, indent=4)
        f.write(contents)
    
def get_bcl_convert_vars(env_vars, sample_sheet, flow_cell):
    vars = env_vars.copy()
    vars['run_name'] = flow_cell
    vars['instrument_platform'] = sample_sheet['instrument_platform'][0]
    vars['instrument_type'] = sample_sheet['instrument_type'][0]
    vars['read1_cycles'] = int(sample_sheet['read1_cycles'][0])
    vars['read2_cycles'] = int(sample_sheet['read2_cycles'][0])
    vars['index1_cycles'] = int(sample_sheet['index1_cycles'][0])
    vars['index2_cycles'] = int(sample_sheet['index2_cycles'][0])
    vars['create_fastq_for_index_reads'] = get_boolean_val(sample_sheet['create_fastq_for_index_reads'][0])
    vars['trim_umi'] = get_boolean_val(sample_sheet['trim_umi'][0])
    vars['override_cycles'] = sample_sheet['override_cycles'][0]
    return vars

def add_lane_to_fastq(file_name):
    if not re.match("^.*L\d{3}.*$", file_name):
        split = re.split("(S\d+)", file_name)
        return ''.join(split[:2]) + '_L001' + split[2]
    return file_name

def get_boolean_val(val):
    match str(val).lower():
        case '1.0' | '1':
            return 1
        case '0.0' | '0':
            return 0 
        case 'true':
            return 'true'
        case 'false':
            return 'false'
        case '' | 'nan':
            return ''
        case _:
            raise ValueError