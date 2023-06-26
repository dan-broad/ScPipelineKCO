version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/kco:bulk_rna_seq/versions/18/plain-WDL/descriptor" as brs 

struct InputArgs {
  String input_bcl_directory
  File bcl_sample_sheet
  File samples_described_file
  String reference_transcriptome
  String organism 
  String output_directory
  Int cycle_number
  String read_pairs_file
}

workflow bulk_rna_seq_pipeline {

  input {
    File sample_sheet

    # BCL CONVERT SETTINGS
    String bcl_convert_version = "3.10.5"
    Boolean no_lane_splitting = true
    Boolean strict_mode = true
    Int minimum_trimmed_read_length = 10
    Int mask_short_reads = 10

    # DE ANALYSIS SETTINGS
    String DE_method
    Float logCPM_threshold = 1
    Float percentage_aligned_threshold = 30 # For sleuth, only consider samples with alignment percentage above this threshold 
    String ensemble_mirror = "useast.ensembl.org"

    # COMPUTE
    Int preemptible = 2
    Int disk_space = 500
    Int num_cpu = 32
    String memory = "120G"
    String zones = "us-east1-d us-west1-a us-west1-b"

  }

  call split_by_flowcell {
    input:
      sample_sheet = sample_sheet
  }

  scatter (flowcell_sample_sheet in split_by_flowcell.sample_sheets) {

    call generate_config {
      input:
        sample_sheet = flowcell_sample_sheet,
        minimum_trimmed_read_length = minimum_trimmed_read_length,
        mask_short_reads = mask_short_reads,
    }

    Map[String, String] config = read_json(generate_config.config)
    call brs.bulk_rna_seq {
      input:
        output_directory = config['output_directory'],
        bcl_sample_sheet = generate_config.bcl_sample_sheet,
        input_data_type = 'bcl',
        input_bcl_directory = config['input_bcl_directory'],
        bcl_convert_version = bcl_convert_version,
        no_lane_splitting = no_lane_splitting,
        strict_mode = strict_mode,
        read_pairs_file = config['read_pairs_file'],

        # KALLISTO
        reference_transcriptome = config['reference_transcriptome'],
        cycle_number = config['cycle_number'],

        # DE
        DE_method = DE_method,
        organism = config['organism'],
        samples_described_file = generate_config.samples_described,
        logCPM_threshold = logCPM_threshold,
        percentage_aligned_threshold = percentage_aligned_threshold,
        ensemble_mirror = ensemble_mirror,

        # COMPUTE
        preemptible = preemptible,
        disk_space = disk_space,
        num_cpu = num_cpu,
        memory = memory,
        zones = zones
    }
  }

}


task split_by_flowcell {
  input {
    File sample_sheet
  }

  command <<<
    set -e
    export TMPDIR=/tmp

    python <<CODE
    import json
    import pandas as pd

    # Load the CSV file
    df = pd.read_csv('~{sample_sheet}', header = 0, dtype = str, index_col = False)
    for seq_dir, group_df in df.groupby('seq_dir'):
      if (group_df['run_pipeline'] == 'TRUE').sum() > 0:
        group_df.to_csv("{}_sample_sheet.csv".format(str(abs(hash(seq_dir)))), index = False)

    CODE
  >>>

  output {
    Array[File] sample_sheets = glob('*_sample_sheet.csv')
  }

  runtime {
    preemptible: 2
    bootDiskSizeGb: 12
    disks: "local-disk 20 HDD"
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    cpu: 16
    zones: "us-east1-d us-west1-a us-west1-b"
    memory: "16G"
  }
}

task generate_config {
  input {
    File sample_sheet
    Int minimum_trimmed_read_length = 10
    Int mask_short_reads = 10
  }

  command <<<
    set -e
    export TMPDIR=/tmp

    python <<CODE
    import json
    import pandas as pd

    # Load the CSV file
    group_df = pd.read_csv('~{sample_sheet}', header = 0, dtype = str, index_col = False)

    # Assert that reference, organism, out_dir, and cycle_number are the same within the group
    assert group_df['reference_transcriptome'].nunique() == 1, "Reference is not the same for seq_dir: {}".format(seq_dir)
    assert group_df['organism'].nunique() == 1, "Organism is not the same for seq_dir: {}".format(seq_dir)
    assert group_df['out_dir'].nunique() == 1, "Out_dir is not the same for seq_dir: {}".format(seq_dir)
    assert group_df['cycle_number'].nunique() == 1, "Cycle_number is not the same for seq_dir: {}".format(seq_dir)

    # Create a new DataFrame with the required columns (sampleid, index, index2)
    index_values = group_df[group_df['run_pipeline'] == 'TRUE'][['sampleid', 'index', 'index2']].copy()
    assert len(index_values) > 0, "At least one sample should have run_pipeline = TRUE"

    # Build bcl convert samplesheet
    pd.DataFrame([
        ['[Header]', '',''],
        ['FileFormatVersion','2',''],
        ['RunName','MyRun', ''],
        ['','',''],
        ['[BCLConvert_Settings]','',''],
        ['SoftwareVersion','3.10.5',''],
        ['MinimumTrimmedReadLength','~{minimum_trimmed_read_length}',''],
        ['MaskShortReads','~{mask_short_reads}',''],
        ['','',''],
        ['[BCLConvert_Data]','',''],
        ['Sample_ID','index','index2']
      ] + list(index_values.values)
    ).to_csv("bcl_sample_sheet.csv", index=False, header=False)

    # Create samples described file
    group_df[group_df['run_pipeline'] == 'TRUE'][['sample_group', 'sampleid']].to_csv("samples_described.tsv", sep='\t', index=None, header=None)

    # Read pairs file location
    out_dir = group_df['out_dir'].iloc[0]
    out_dir = out_dir + '/' if out_dir[-1] != '/' else out_dir
    read_pairs_file = "{}{}_read_pairs.tsv".format(out_dir, str(abs(hash(group_df['seq_dir'].iloc[0]))))

    # Create a dictionary entry for the group
    config_dict = {
      'input_bcl_directory': group_df['seq_dir'].iloc[0],
      'reference_transcriptome': group_df['reference_transcriptome'].iloc[0],
      'organism': group_df['organism'].iloc[0],
      'output_directory': group_df['out_dir'].iloc[0],
      'cycle_number': group_df['cycle_number'].iloc[0],
      'read_pairs_file': read_pairs_file
    }
    config_name = "config.json"
    with open(config_name, "w") as outfile:
      json.dump(config_dict, outfile)
    CODE
  
  >>>

  output {
    File config = 'config.json'
    File bcl_sample_sheet = 'bcl_sample_sheet.csv'
    File samples_described = 'samples_described.tsv'
  }

  
  runtime {
    preemptible: 2
    bootDiskSizeGb: 12
    disks: "local-disk 20 HDD"
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    cpu: 16
    zones: "us-east1-d us-west1-a us-west1-b"
    memory: "16G"
  }
}
