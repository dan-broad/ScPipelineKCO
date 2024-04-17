import logging
import os
import re
import time
import sys
import subprocess
import threading
import json
from string import Template
import firecloud.api as fapi
import pandas as pd
from pathlib import Path
from consts import TERRA_POLL_SPACER, TERRA_TIMEOUT, RNA, ATAC, GEX_I7_INDEX_KEY

alto_lock = threading.Lock()


def build_directories(basedir):
    directories = {
        'fastqs': basedir + "/fastqs",
        'counts': basedir + "/counts",
        'results': basedir + "/cumulus",
        'cellranger_arc': basedir + "/cellranger_arc",
        'cellbender': basedir + "/cellbender",
        'cellbender_results': basedir + "/cellbender_cumulus",
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
        'cellbender': gcp_basedir + "/cellbender_" + project,
        'cellbender_results': gcp_basedir + "/cellbender_cumulus_" + project,
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
    insert_cellbender_defaults(sample_tracking)
    
    for _, row in sample_tracking.iterrows():
        learning_rate = "%f" % row['cellbender_learning_rate'] # gets rid of scientific notation for floats, temp fix
        
        sample_dict[row['sampleid']].append(row['Sample'])
        mkfastq_dict[row['Sample']] = [row['Lane'], row['Index'], row['reference'], row['chemistry'], row['method']]
        cumulus_dict[row['sampleid']] = [row['min_umis'], row['min_genes'], row['percent_mito']]
        cellbender_dict[row['sampleid']] = [row['cellbender_expected_cells'],
                                            row['cellbender_total_droplets_included'],
                                            learning_rate,
                                            row['cellbender_force_cell_umi_prior'],
                                            row['cellbender_force_empty_umi_prior']]
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
    file_dir = os.path.dirname(os.path.realpath(__file__))
    with open(f'{file_dir}/templates/bcl_convert_gex_sample_sheet_template.csv') as f:
        template = Template(f.read())
        sample_sheet = template.safe_substitute(env_vars)

    columns = ["sampleid", "Index"]
    columns = columns + ["Lane"] if not env_vars.get("no_lane_splitting") else columns
    samples_with_indices = sample_tracking[columns]

    if not env_vars.get("no_lane_splitting"):
        samples_with_indices = apply_lane_splits(samples_with_indices, num_lanes)
        samples_with_indices = samples_with_indices[['Lane', 'sampleid', 'Index']] #reorder columns

    samples_with_indices.rename(columns={"sampleid": "Sample_ID", "Index": "index"}, inplace=True)
    replace_index(samples_with_indices, RNA, env_vars["gex_i5_index_key"])
    return sample_sheet + samples_with_indices.to_csv(index=False,)

def get_atac_sample_sheet(env_vars, sample_tracking, num_lanes):
    file_dir = os.path.dirname(os.path.realpath(__file__))
    with open(f'{file_dir}/templates/bcl_convert_atac_sample_sheet_template.csv') as f:
        template = Template(f.read())
        sample_sheet = template.safe_substitute(env_vars)

    columns = ["Lane", "sampleid", "Index"]
    samples_with_indices = sample_tracking.get(columns)
    samples_with_indices.rename(columns={"sampleid": "Sample_ID", "Index": "index"}, inplace=True)
    replace_index(samples_with_indices, ATAC)
    flattend_sample_indices = samples_with_indices.explode("index")

    if env_vars.get("no_lane_splitting"):
        flattend_sample_indices.drop(columns='Lane', inplace=True)
    else:
        flattend_sample_indices = apply_lane_splits(flattend_sample_indices, num_lanes)
        flattend_sample_indices = flattend_sample_indices.get(['Lane', 'Sample_ID', 'index', 'index2']) #reorder columns

    return sample_sheet + flattend_sample_indices.to_csv(index=False)

def get_library_indices():
    file_dir = os.path.dirname(os.path.realpath(__file__))
    with open(f"{file_dir}/index_kits/Dual_Index_Kit_TT_Set_A.json") as gex_indices, \
        open(f"{file_dir}/index_kits/Single_Index_Kit_N_Set_A_Reformatted.json") as atac_indices:
        return {
            RNA: json.load(gex_indices),
            ATAC: json.load(atac_indices)
        }

def replace_index(samples, method, i5_index_key=None):
    indices = get_library_indices()[method]
    if method == RNA:
        samples.insert(len(samples.columns), 'index2', pd.NA)

    for i, r in samples.iterrows():
        index_oligonucleotide = indices[r['index']]
        if method == RNA:
            samples.loc[i, 'index'] = index_oligonucleotide[GEX_I7_INDEX_KEY]
            samples.loc[i, 'index2'] = index_oligonucleotide[i5_index_key]
        elif method == ATAC:
            samples.loc[i, 'index'] = index_oligonucleotide

def apply_lane_splits(sample_tracking, num_lanes):
    for _, sample in sample_tracking.iterrows():
        lane_val = str(sample['Lane'])

        if len(lane_val) == 1 and lane_val.isdecimal(): # no need to split
            sample['Lane'] = sample['Lane'].to_list()
        elif '-' in lane_val:
            start, end = lane_val.split('-')
            sample['Lane'] = list(range(int(start), int(end)+1))
        elif lane_val == '*':
            sample['Lane'] = list(range(1, num_lanes+1))

    return sample_tracking.explode('Lane')

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
    bcl_convert_vars = env_vars.copy()
    override_cycles = sample_sheet['override_cycles'][0]
    read_cycles = re.findall("Y\d+", override_cycles)
    index_cycles = re.findall("[I|U]\d+", override_cycles)

    bcl_convert_vars['run_name'] = flow_cell
    bcl_convert_vars['instrument_platform'] = sample_sheet['instrument_platform'][0]
    bcl_convert_vars['instrument_type'] = sample_sheet['instrument_type'][0]
    bcl_convert_vars['read1_cycles'] = int(read_cycles[0][1:])
    bcl_convert_vars['read2_cycles'] = int(read_cycles[1][1:])
    bcl_convert_vars['index1_cycles'] = int(index_cycles[0][1:])
    bcl_convert_vars['index2_cycles'] = int(index_cycles[1][1:])
    bcl_convert_vars['create_fastq_for_index_reads'] = get_boolean_val(sample_sheet['create_fastq_for_index_reads'][0])
    bcl_convert_vars['trim_umi'] = get_boolean_val(sample_sheet['trim_umi'][0])
    bcl_convert_vars['override_cycles'] = override_cycles
    return bcl_convert_vars

def add_lane_to_fastq(file_name):
    if not re.match("^.*L\d{3}.*$", file_name):
        split = re.split("(S\d+)", file_name)
        return ''.join(split[:2]) + '_L001' + split[2]
    return file_name

# def get_boolean_val(val):
#     match str(val).lower():
#         case '1.0' | '1':
#             return 1
#         case '0.0' | '0':
#             return 0 
#         case 'true':
#             return 'true'
#         case 'false':
#             return 'false'
#         case '' | 'nan':
#             return ''
#         case _:
#             raise ValueError

def get_boolean_val(val):
    val = str(val).lower():
    if val in ('1.0', '1'):
        return 1
    elif val in ('0.0', '0'):
        return 0 
    elif val == 'true':
        return 'true'
    elif val == 'false':
        return 'false'
    elif val in ('', 'nan'):
        return ''
    else:
        raise ValueError

def get_cellbender_inputs_template(version):
    parent_dir = Path(__file__).parent.resolve()

    if version.startswith("0.3"):
        with open(f'{parent_dir}/templates/cellbender_v3_input_template.json') as f:
            template = f.read()
    else:
        with open(f'{parent_dir}/templates/cellbender_v2_input_template.json') as f:
            template = f.read()
    return template

def insert_cellbender_defaults(sample_sheet):
    cellbender_defaults = {
        "cellbender_learning_rate": 0.00005,
    }

    sample_sheet.fillna(value=cellbender_defaults, inplace=True)
