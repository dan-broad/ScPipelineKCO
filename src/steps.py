import logging
import os
from datetime import datetime
from pathlib import Path
from utils import *
from consts import MULTIOME

def upload_bcl_convert_input(sample_tracking, buckets, directories, env_vars):
    logging.info("STEP 1 | Upload bcl convert samplesheet to Google Cloud Storage Bucket. ")

    fc = sample_tracking['flowcell'].to_list()[0]
    samples_in_dir = sample_tracking.where(sample_tracking['flowcell'] == fc).dropna(how='all')
    method = samples_in_dir['sub_method'].to_list()[0]
    directory = f"{directories['bcl_convert']}/{fc}/{method}"
    bucket_path = f'{buckets["bcl_convert"]}/inputs/{fc}'
    os.makedirs(directory, exist_ok=True)
    input_paths = {"bucket": bucket_path, "directory": directory}

    bcl_convert_vars = get_bcl_convert_vars(env_vars, samples_in_dir, fc)
    file_name = "sample_sheet.csv"
    path = f"{directory}/{file_name}"
    create_bcl_convert_sample_sheet(path, method, bcl_convert_vars, samples_in_dir)
    input_dir = samples_in_dir['seq_dir'][0]
    create_bcl_convert_params(f"{directory}/inputs.json", env_vars, input_dir, f"{buckets['fastqs']}/{method}", f"{bucket_path}/{file_name}")

    upload_bcl_convert_file = "%s/upload_%s.sh" % (directories['bcl_convert'], datetime.now().isoformat())
    with open(upload_bcl_convert_file, "w") as f:
        f.write(f"gsutil cp {input_paths['directory']}/* {input_paths['bucket']}\n")

    bash_execute_file(upload_bcl_convert_file)
    return input_paths

def run_bcl_convert(directories, buckets, sample_paths, bcl_convert_method, bcl_convert_workspace):
    run_alto_file = f"{directories['bcl_convert']}/run_bcl_convert_workflow_{datetime.now().isoformat()}.sh"

    with open(run_alto_file, "w") as f:
        directory_path = Path(sample_paths['directory'])
        bucket_fastq_path = Path(buckets["fastqs"])
        method = directory_path.parent.name
        input_file = f"{directory_path}/inputs.json"
        f.write(f"alto terra run -m {bcl_convert_method} -i {input_file} -w {bcl_convert_workspace} --bucket-folder {bucket_fastq_path.name}/{method} --no-cache\n")

    logging.info("STEP 2 | Initiate Terra bcl_convert pipeline via alto. ")
    execute_alto_command(run_alto_file)

def get_fastq_paths(directories, buckets, sample_tracking):
    get_fastq_paths_file = f"{directories['bcl_convert']}/get_fastqs_{datetime.now().isoformat()}.sh"
    fastq_paths_file = f"{directories['fastqs']}/fastq_dirs_{datetime.now().isoformat()}.txt" # using runtime here to not create multiple files that we have to write then read 
    with open(get_fastq_paths_file, "w") as f:
        for _, r in sample_tracking.iterrows():
            run_id = os.path.basename(r['seq_dir'])
            sample_id = r['sampleid']
            base_path = f"{buckets['fastqs']}/{r['sub_method']}/{run_id}_fastqs/sample_fastqs"
            cmd = f"gsutil ls '{base_path}/{sample_id}**' >> {fastq_paths_file}"
            f.write(f"{cmd}\n")
    bash_execute_file(get_fastq_paths_file)
    return fastq_paths_file

def move_fastqs_to_sample_dir(directories, buckets, sample_tracking):
    fastqs_paths_file = open(get_fastq_paths(directories, buckets, sample_tracking))
    fastq_paths = fastqs_paths_file.readlines()
    fastqs_paths_file.close()
    
    move_fastqs_file = f"{directories['bcl_convert']}/move_fastqs_{datetime.now().isoformat()}.sh"
    with open(move_fastqs_file, "w") as f:
        for fp in fastq_paths:
            fp = fp.strip()
            p = Path(fp)
            base_path = str(p.parent.parent)
            base_path = base_path.replace("gs:/", "gs://") if "gs://" not in base_path else base_path
            fastq_name = p.name
            fastq_name_with_lane = add_lane_to_fastq(fastq_name)
            sample_id = fastq_name.split("_")[0]

            #update fastq sample ids to match ids expected by cellranger arc
            if sample_tracking[sample_tracking['sampleid'] == sample_id]['method'].to_list()[0] == MULTIOME:
                submethod = p.parent.parent.parent.parent.name
                fastq_name_with_lane = f"{sample_id}_{submethod}_{'_'.join(fastq_name_with_lane.split('_')[1:])}"

            f.write(f"gsutil mv {fp} {base_path}/{sample_id}/{fastq_name_with_lane}\n")

    logging.info("Moving FASTQs to sample directory. ")
    bash_execute_file(move_fastqs_file)

@DeprecationWarning
def upload_cellranger_mkfastq_input(buckets, directories, sample_tracking, cellranger_version, cellranger_atac_version, diskspace, memory):
    run_id = os.path.basename(sample_tracking['seq_dir'].tolist()[0])
    fastq_flowcell_bucket = "%s/%s_fastqs" % (buckets['fastqs'], run_id)
    fastq_flowcell_dir = "%s/%s" % (directories['fastqs'], run_id)
    os.makedirs(fastq_flowcell_dir, exist_ok=True)

    samplesheet_file = '%s/samplesheet_cellranger_mkfastq.csv' % fastq_flowcell_dir
    samplesheet = sample_tracking[['Lane', 'Sample', 'Index', 'reference', 'chemistry', 'seq_dir', 'method']]
    samplesheet.to_csv(samplesheet_file, index=False,
                       header=['Lane', 'Sample', 'Index', 'Reference', 'Chemistry', 'Flowcell', 'DataType'])

    input_cellranger_file = "%s/input_cellranger_mkfastq.json" % fastq_flowcell_dir

    with open(input_cellranger_file, "w") as f:
        f.write("{\n")
        f.write("\t\"cellranger_workflow.input_csv_file\" : \"%s/samplesheet_cellranger_mkfastq.csv\",\n" %
                fastq_flowcell_bucket)
        f.write("\t\"cellranger_workflow.output_directory\" : \"%s\",\n" % buckets['fastqs'])
        f.write("\t\"cellranger_workflow.cellranger_version\" : \"%s\",\n" % cellranger_version)
        f.write("\t\"cellranger_workflow.cellranger_atac_version\" : \"%s\",\n" % cellranger_atac_version)
        f.write("\t\"cellranger_workflow.run_mkfastq\" : true,\n")
        f.write("\t\"cellranger_workflow.run_count\" : false,\n")
        f.write("\t\"cellranger_workflow.mkfastq_disk_space\" : %s,\n" % diskspace)
        f.write("\t\"cellranger_workflow.memory\" : \"%s\",\n" % memory)
        f.write("\t\"cellranger_workflow.mkfastq_docker_registry\" : \"gcr.io/microbiome-xavier\"\n")
        f.write("}\n")

    logging.info("STEP 1 | Upload cellranger samplesheet and input file to Google Cloud Storage Bucket. ")
    upload_cellranger_file = "%s/upload_cellranger_mkfastq.sh" % fastq_flowcell_dir
    with open(upload_cellranger_file, "w") as f:
        f.write("gsutil cp %s %s/\n" % (samplesheet_file, fastq_flowcell_bucket))
        f.write("gsutil cp %s %s/\n" % (input_cellranger_file, fastq_flowcell_bucket))
    bash_execute_file(upload_cellranger_file)

@DeprecationWarning
def run_cellranger_mkfastq(directories, sample_tracking, alto_workspace, alto_method, alto_fastqs_folder):
    run_id = os.path.basename(sample_tracking['seq_dir'].iloc[0])
    alto_flowcell_bucket = "%s/%s" % (alto_fastqs_folder, run_id)
    fastq_flowcell_dir = directories['fastqs'] + "/%s" % run_id

    run_alto_file = "%s/run_alto_cellranger_workflow.sh" % fastq_flowcell_dir

    open(run_alto_file, "w").close()
    bash_alto = open(run_alto_file, "a")
    input_cellranger_file = "%s/input_cellranger_mkfastq.json" % fastq_flowcell_dir
    bash_alto.write("alto terra run -m %s -i %s -w %s --bucket-folder %s --no-cache\n" % (
        alto_method, input_cellranger_file, alto_workspace, alto_flowcell_bucket))
    bash_alto.close()

    logging.info("STEP 2 | Initiate Terra cellranger_workflow pipeline via alto. ")
    execute_alto_command(run_alto_file)


def upload_cellranger_count_input(buckets, directories, sample_dicts, sample_tracking, cellranger_version, cellranger_atac_version):
    fastqs_bucket = buckets['fastqs']
    counts_bucket = buckets['counts']
    counts_dir = directories['counts']
    sample_dict = sample_dicts['sample']
    mkfastq_dict = sample_dicts['mkfastq']
    cellranger_dict = sample_dicts['cellranger']

    run_id = os.path.basename(sample_tracking['seq_dir'].tolist()[0])
    flowcell = sample_tracking['flowcell'].iloc[0]
    method = sample_tracking['sub_method'].iloc[0]
    fastq_flowcell_bucket = "%s/%s/%s_fastqs/sample_fastqs" % (fastqs_bucket, method, run_id)

    for sample_id in sample_dict.keys():
        os.makedirs("%s/%s" % (counts_dir, sample_id), exist_ok=True)
        samplesheet_cellranger_file = "%s/%s/samplesheet_cellranger.csv" % (counts_dir, sample_id)

        with open(samplesheet_cellranger_file, "w") as f:
            f.write("Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType\n")
            for sample in sample_dict[sample_id]:
                fastq_sample_bucket = "%s/%s/" % (fastq_flowcell_bucket, sample_id)
                f.write("%s,%s,%s,%s,%s,%s,%s\n" % (sample, mkfastq_dict[sample][2], fastq_sample_bucket,
                                                 mkfastq_dict[sample][0], mkfastq_dict[sample][1],
                                                 mkfastq_dict[sample][3],mkfastq_dict[sample][4]))

    for sample_id in sample_dict.keys():
        input_cellranger_file = "%s/%s/input_cellranger.json" % (counts_dir, sample_id)

        with open(input_cellranger_file, "w") as f:
            f.write("{\n")
            f.write("\t\"cellranger_workflow.input_csv_file\" : \"%s/%s/samplesheet_cellranger.csv\",\n" % (
                counts_bucket, sample_id))
            f.write("\t\"cellranger_workflow.output_directory\" : \"%s\",\n" % counts_bucket)
            f.write("\t\"cellranger_workflow.cellranger_version\" : \"%s\",\n" % cellranger_version)
            f.write("\t\"cellranger_workflow.cellranger_atac_version\" : \"%s\",\n" % cellranger_atac_version)
            f.write("\t\"cellranger_workflow.run_mkfastq\" : false,\n")
            f.write("\t\"cellranger_workflow.run_count\" : true,\n")
            f.write("\t\"cellranger_workflow.mkfastq_docker_registry\" : \"gcr.io/microbiome-xavier\",\n")
            f.write("\t\"cellranger_workflow.include_introns\" : %s\n" % str(cellranger_dict[sample_id][0]).lower())
            f.write("}\n")

    logging.info("STEP 3 | Upload cellranger samplesheet and input file to Google Cloud Storage Bucket. ")
    uploadcellranger_file = "%s/uploadcellranger_%s.sh" % (counts_dir, flowcell)
    with open(uploadcellranger_file, "w") as f:
        for sample_id in sample_dict.keys():
            samplesheet_cellranger_file = "%s/%s/samplesheet_cellranger.csv" % (counts_dir, sample_id)
            input_cellranger_file = "%s/%s/input_cellranger.json" % (counts_dir, sample_id)
            f.write("gsutil cp %s %s/%s/\n" % (samplesheet_cellranger_file, counts_bucket, sample_id))
            f.write("gsutil cp %s %s/%s/\n" % (input_cellranger_file, counts_bucket, sample_id))
    bash_execute_file(uploadcellranger_file)


def run_cellranger_count(directories, sample_dicts, sample_tracking, alto_workspace, alto_method, alto_counts_folder):
    sample_dict = sample_dicts['sample']
    counts_dir = directories['counts']

    run_alto_file = "%s/run_alto_cellranger_workflow_%s.sh" % (counts_dir, sample_tracking['flowcell'].iloc[0])

    open(run_alto_file, "w").close()
    bash_alto = open(run_alto_file, "a")
    for sampleid in sample_dict.keys():
        input_cellranger_file = "%s/%s/input_cellranger.json" % (counts_dir, sampleid)
        bash_alto.write("alto terra run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (
            alto_method, input_cellranger_file, alto_workspace, alto_counts_folder, sampleid))
    bash_alto.close()

    logging.info("STEP 4 | Initiate Terra cellranger_workflow pipeline via alto. ")
    execute_alto_command(run_alto_file)


def upload_cumulus_samplesheet(buckets, directories, sample_dicts, sample_tracking, count_matrix_name):
    sample_dict = sample_dicts['sample']
    cumulusdict = sample_dicts['cumulus']
    results_dir = directories['results']
    resultsbucket = buckets['results']
    parent_dir = Path(__file__).parent.resolve()
    
    flowcell = sample_tracking['flowcell'].iloc[0]

    for sampleid in sample_dict.keys():
        if not os.path.isdir("%s/%s" % (results_dir, sampleid)):
            os.mkdir("%s/%s" % (results_dir, sampleid))
        samplesheet_cumulus_file = "%s/%s/samplesheet_cumulus.csv" % (results_dir, sampleid)

        with open(samplesheet_cumulus_file, "w") as f:
            f.write("Sample,Location\n")
            for sample in sample_dict[sampleid]:
                f.write("%s,%s/%s/%s\n" % (sample, buckets['counts'], sample, count_matrix_name))

    
    for sampleid in sample_dict.keys():
        input_cumulus_file = "%s/%s/input_cumulus.json" % (results_dir, sampleid)
        samplesheet_cumulus = "%s/%s/samplesheet_cumulus.csv" % (resultsbucket, sampleid)

        with open(f'{parent_dir}/templates/cumulus_input_template.json') as f:
            template = f.read()
        template = template.replace('{input_file}', samplesheet_cumulus) \
            .replace('{output_directory}', "%s/" % resultsbucket) \
            .replace('{output_name}', sampleid) \
            .replace('"{min_umis}"', str(cumulusdict[sampleid][0])) \
            .replace('"{min_genes}"', str(cumulusdict[sampleid][1])) \
            .replace('"{percent_mito}"', str(cumulusdict[sampleid][2]))

        with open(input_cumulus_file, "w") as f:
            f.write(template)

    logging.info("STEP 5 | Upload cumulus samplesheet and input file to Google Cloud Storage Bucket. ")
    upload_cumulus_file = "%s/uploadcumulus_%s.sh" % (results_dir, flowcell)
    with open(upload_cumulus_file, "w") as f:
        for sampleid in sample_dict.keys():
            samplesheet_cumulus_file = "%s/%s/samplesheet_cumulus.csv" % (results_dir, sampleid)
            input_cumulus_file = "%s/%s/input_cumulus.json" % (results_dir, sampleid)
            f.write("gsutil cp %s %s/%s/\n" % (samplesheet_cumulus_file, resultsbucket, sampleid))
            f.write("gsutil cp %s %s/%s/\n" % (input_cumulus_file, resultsbucket, sampleid))
    bash_execute_file(upload_cumulus_file)


def run_cumulus(directories, sample_dicts, sample_tracking, alto_workspace, alto_method, alto_results_folder):
    sample_dict = sample_dicts['sample']
    results_dir = directories['results']

    run_alto_file = "%s/run_alto_cumulus_%s.sh" % (results_dir, sample_tracking['flowcell'].iloc[0])

    open(run_alto_file, "w").close()
    bash_alto = open(run_alto_file, "a")
    for sampleid in sample_dict.keys():
        input_cumulus_file = "%s/%s/input_cumulus.json" % (results_dir, sampleid)
        bash_alto.write("alto terra run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (
            alto_method, input_cumulus_file, alto_workspace, alto_results_folder, sampleid))
    bash_alto.close()

    logging.info("STEP 6 | Initiate Terra cumulus pipeline via alto. ")
    execute_alto_command(run_alto_file)


def upload_cell_bender_input(buckets, directories, sample_dicts, sample_tracking, count_matrix_name, version):
    sample_dict = sample_dicts['sample']
    cellbenderdict = sample_dicts['cellbender']
    cellbender_dir = directories['cellbender']
    cellbenderbucket = buckets['cellbender']
    countsbucket = buckets['counts']

    for sampleid in sample_dict.keys():
        if not os.path.isdir("%s/%s" % (cellbender_dir, sampleid)):
            os.mkdir("%s/%s" % (cellbender_dir, sampleid))

        for sample in sample_dict[sampleid]:
            if not os.path.isdir("%s/%s" % (cellbender_dir, sampleid)):
                os.mkdir("%s/%s" % (cellbender_dir, sampleid))
            input_cellbender_file = "%s/%s/input_cellbender.json" % (cellbender_dir, sampleid)
            template = get_cellbender_inputs_template(version)

            template = template.replace('{sample_name}', str(sampleid)) \
                .replace('{input_dir}', "%s/%s/%s" % (countsbucket, sampleid, count_matrix_name)) \
                .replace('{output_dir}', "%s/%s" % (cellbenderbucket, sampleid)) \
                .replace('{cellbender_version}', version) \
                .replace('{learning_rate}', np.format_float_positional(cellbenderdict[sampleid][2]))

            optional_params = {'"{total_droplets_included}"': str(cellbenderdict[sampleid][1]),
                               '"{expected_cells}"': str(cellbenderdict[sampleid][0]),
                               '"{force_cell_umi_prior}"': str(cellbenderdict[sampleid][3]),
                               '"{force_empty_umi_prior}"': str(cellbenderdict[sampleid][4])
                            }
            
            for replace_from, replace_to in optional_params.items():
                if pd.notna(replace_to) and replace_to.lower() not in ['nan', '']:
                    template = template.replace(replace_from, replace_to)
                else:
                    template = template.replace(replace_from, 'null')

            template = json.loads(template)
            params = {}
            for k, v in template.items():
                if not v is None:
                    params[k] = v
                if 'learning_rate' in k:
                    params[k] = np.format_float_positional(v) # get float out of scientific notation

            with open(input_cellbender_file, "w") as f:
                json.dump(params, f)

    logging.info("STEP 7 | Upload cellbender input file to Google Cloud Storage Bucket. ")
    uploadcellbender_file = "%s/uploadcellbender_%s.sh" % (cellbender_dir, sample_tracking['flowcell'].iloc[0])
    with open(uploadcellbender_file, "w") as f:
        for sampleid in sample_dict.keys():
            for sample in sample_dict[sampleid]:
                input_cellbender_file = "%s/%s/input_cellbender.json" % (cellbender_dir, sampleid)
                f.write("gsutil cp %s %s/%s/\n" % (input_cellbender_file, cellbenderbucket, sampleid))
    bash_execute_file(uploadcellbender_file)


def run_cellbender(directories, sample_dicts, sample_tracking, alto_workspace, alto_method, alto_cellbender_folder):
    sample_dict = sample_dicts['sample']
    cellbender_dir = directories['cellbender']

    run_alto_file = "%s/run_alto_cellbender_%s.sh" % (cellbender_dir, sample_tracking['flowcell'].iloc[0])

    open(run_alto_file, "w").close()
    bash_alto = open(run_alto_file, "a")
    for sampleid in sample_dict.keys():
        for sample in sample_dict[sampleid]:
            input_cellbender_file = "%s/%s/input_cellbender.json" % (cellbender_dir, sampleid)
            bash_alto.write("alto terra run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (
                alto_method, input_cellbender_file, alto_workspace, alto_cellbender_folder, sampleid))
    bash_alto.close()

    logging.info("STEP 8 | Initiate Terra remove-background pipeline via alto. ")
    execute_alto_command(run_alto_file)


def upload_post_cellbender_cumulus_input(buckets, directories, sample_dicts, sample_tracking, cellbender_matrix_name):
    sample_dict = sample_dicts['sample']
    cumulusdict = sample_dicts['cumulus']
    cellbenderbucket = buckets['cellbender']
    cellbender_resultsbucket = buckets['cellbender_results']
    cellbender_results_dir = directories['cellbender_results']
    parent_dir = Path(__file__).parent.resolve()

    for sampleid in sample_dict.keys():
        if not os.path.isdir("%s/%s" % (cellbender_results_dir, sampleid)):
            os.mkdir("%s/%s" % (cellbender_results_dir, sampleid))
        samplesheet_cellbender_cumulus_file = "%s/%s/samplesheet_cellbender_cumulus.csv" % (
            cellbender_results_dir, sampleid)

        with open(samplesheet_cellbender_cumulus_file, "w") as f:
            f.write("Sample,Location\n")
            for sample in sample_dict[sampleid]:
                f.write("%s,%s/%s/%s/%s_%s\n" % (
                    sample, cellbenderbucket, sampleid, sampleid, sampleid, cellbender_matrix_name))

    for sampleid in sample_dict.keys():
        input_cellbender_cumulus_file = "%s/%s/input_cumulus.json" % (cellbender_results_dir, sampleid)
        samplesheet_cellbender = "%s/%s/samplesheet_cellbender_cumulus.csv" % (cellbender_resultsbucket, sampleid)

        with open(f'{parent_dir}/templates/cumulus_input_template.json') as f:
            template = f.read()
        template = template.replace('{input_file}', samplesheet_cellbender) \
            .replace('{output_directory}', "%s/" % cellbender_resultsbucket) \
            .replace('{output_name}', sampleid) \
            .replace('"{min_umis}"', str(cumulusdict[sampleid][0])) \
            .replace('"{min_genes}"', str(cumulusdict[sampleid][1])) \
            .replace('"{percent_mito}"', str(cumulusdict[sampleid][2])) \

        with open(input_cellbender_cumulus_file, "w") as f:
            f.write(template)

    logging.info("STEP 9 | Upload post-cellbender cumulus samplesheet and input file to Google Cloud Storage Bucket. ")
    uploadcellbendercumulus_file = "%s/uploadcellbendercumulus_%s.sh" % (
        cellbender_results_dir, sample_tracking['flowcell'].iloc[0])
    with open(uploadcellbendercumulus_file, "w") as f:
        for sampleid in sample_dict.keys():
            samplesheet_cellbender_cumulus_file = "%s/%s/samplesheet_cellbender_cumulus.csv" % (
                cellbender_results_dir, sampleid)
            input_cellbender_cumulus_file = "%s/%s/input_cumulus.json" % (cellbender_results_dir, sampleid)
            f.write("gsutil cp %s %s/%s/\n" % (samplesheet_cellbender_cumulus_file, cellbender_resultsbucket, sampleid))
            f.write("gsutil cp %s %s/%s/\n" % (input_cellbender_cumulus_file, cellbender_resultsbucket, sampleid))
    bash_execute_file(uploadcellbendercumulus_file)


def run_cumulus_post_cellbender(directories, sample_dicts, sample_tracking, alto_workspace, alto_method, alto_results_folder):
    sample_dict = sample_dicts['sample']
    cellbender_results_dir = directories['cellbender_results']

    run_alto_file = "%s/run_alto_cellbender_cumulus_%s.sh" % (cellbender_results_dir, sample_tracking['flowcell'].iloc[0])

    open(run_alto_file, "w").close()
    bash_alto = open(run_alto_file, "a")
    for sampleid in sample_dict.keys():
        input_cellbender_cumulus_file = "%s/%s/input_cumulus.json" % (cellbender_results_dir, sampleid)
        bash_alto.write("alto terra run -m %s -i %s -w %s --bucket-folder %s/%s --no-cache\n" % (
            alto_method, input_cellbender_cumulus_file, alto_workspace, alto_results_folder, sampleid))
    bash_alto.close()

    logging.info("STEP 10 | Initiate Terra cumulus pipeline via alto. ")
    execute_alto_command(run_alto_file)


def upload_cellranger_arc_samplesheet(buckets, directories, sample_tracking, cellranger_arc_version,
                                      mkfastq_disk_space, mkfastq_memory, steps_to_run):
    arc_dir = directories['cellranger_arc']
    arc_bucket = buckets['cellranger_arc']
    parent_dir = Path(__file__).parent.resolve()

    samplesheet_arc_file = f"{arc_dir}/arc/samplesheet_arc.csv"
    samplesheet_arc_gcp_file = f"{arc_bucket}/arc/samplesheet_arc.csv"

    input_arc_file = f"{arc_dir}/arc/input_arc.json"
    input_arc_gcp_file = f"{arc_bucket}/arc/input_arc.json"

    if not os.path.isdir(f"{arc_dir}/arc"):
        os.mkdir(f"{arc_dir}/arc")

    with open(samplesheet_arc_file, "w") as f:
        f.write("Sample,Reference,Flowcell,Lane,Index,DataType,Link\n")
        for _, sample in sample_tracking.iterrows():
            sample_id = f"{sample['Sample']}_{sample['sub_method']}"
            f.write(f"{sample_id},{sample['reference']},{sample['fastq_dir']},{sample['Lane']},{sample['Index']},{sample['sub_method']},{sample['Sample']}\n")

    include_introns = set(sample_tracking['introns'])
    if len(include_introns) != 1:
        logging.error("Unable to run samples with introns included and without in the same run. Exiting.")
        exit(1)
    include_introns = include_introns.pop()
    run_mkfastq = "MKFASTQ" in steps_to_run
    run_count = "COUNT" in steps_to_run

    
    with open(f'{parent_dir}/templates/cellranger_arc_input_template.json') as f:
        template = f.read().replace('{input_csv}', samplesheet_arc_gcp_file) \
            .replace('{output_dir}', f"{arc_bucket}/output/") \
            .replace('"{include_introns}"', f'{str(include_introns).lower()}') \
            .replace('{cellranger_arc_version}', f'{cellranger_arc_version}') \
            .replace('"{mkfastq_disk_space}"', f'{mkfastq_disk_space}') \
            .replace('{memory}', f'{mkfastq_memory}') \
            .replace('"{run_mkfastq}"', f'{str(run_mkfastq).lower()}') \
            .replace('"{run_count}"', f'{str(run_count).lower()}')

    with open(input_arc_file, "w") as f:
        f.write(template)

    logging.info("MULTIOME | Upload cellranger arc samplesheet and input file to Google Cloud Storage Bucket. ")
    upload_arc_file = f"{arc_dir}/upload_arc.sh"
    with open(upload_arc_file, "w") as f:
        f.write(f"gsutil cp {samplesheet_arc_file} {samplesheet_arc_gcp_file}\n")
        f.write(f"gsutil cp {input_arc_file} {input_arc_gcp_file}\n")

    bash_execute_file(upload_arc_file)


def run_cellranger_arc(buckets, directories, alto_method, alto_workspace):
    arc_dir = directories['cellranger_arc']
    arc_bucket = buckets['cellranger_arc']

    run_alto_file = f"{arc_dir}/run_alto_cellranger_arc.sh"

    with open(run_alto_file, "w") as f:
        input_arc_file = f"{arc_dir}/arc/input_arc.json"
        f.write(f"alto terra run -m {alto_method} -i {input_arc_file} -w {alto_workspace} --bucket-folder {arc_bucket}\n")

    logging.info("STEP 6 | Initiate Terra cumulus pipeline via alto. ")
    execute_alto_command(run_alto_file)
