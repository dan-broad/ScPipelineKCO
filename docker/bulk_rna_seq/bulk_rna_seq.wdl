version 1.0


workflow bulk_rna_seq {

  input {
    String output_directory
    
    # BCL Convert args
    File? bcl_sample_sheet
    String input_data_type
    String input_bcl_directory = ''
    String bcl_convert_version = "3.10.5"
    Boolean no_lane_splitting = true
    Boolean strict_mode = true
    String read_pairs_file = ''

    # KALLISTO
    File reference_transcriptome
    Int cycle_number = 31

    # DE
    String DE_method
    String organism
    String samples_described_file
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

  String output_directory_stripped = sub(output_directory, "/+$", "")

  call run_bcl_convert_and_define_read_pairs {
    input:
      input_data_type = input_data_type,
      bcl_sample_sheet = bcl_sample_sheet,
      input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
      output_directory = output_directory_stripped,
      bcl_convert_version = bcl_convert_version,
      no_lane_splitting = no_lane_splitting,
      strict_mode = strict_mode,
      read_pairs_file = read_pairs_file,
      
      zones = zones,
      num_cpu = num_cpu,
      memory = memory,
      disk_space = disk_space,
      preemptible = preemptible,
  }

  # build reference transcriptome index for kallisto
  call build_kallisto_index {
    input:
      reference_transcriptome = reference_transcriptome,
      cycle_number = cycle_number
  }

  # get pairs of reads
  Array[Pair[String, String]] zipped_read_pairs = read_map(run_bcl_convert_and_define_read_pairs.read_pairs_tsv)

  # for each read pair ...
  scatter (pair in zipped_read_pairs) {
    # define files matching read pair
    File zipped_R1 = "~{pair.left}"
    File zipped_R2 = "~{pair.right}"

    # filter read name to contain only sample information
    String sample_name = sub(sub(basename(zipped_R1), ".fastq.gz", ""), "_R[0-9]", "")

    # unzip each file (if zipped) -- this is necessary for fastqc
    call unzip_file {
      input:
        zipped_file_1 = zipped_R1,
        zipped_file_2 = zipped_R2
    }

    # perform fastqc
    call perform_fastqc {
      input:
        output_directory = output_directory_stripped,
        file_1 = unzip_file.file_1,
        file_2 = unzip_file.file_2,
        sample_name = sample_name
    }

    # get transcript quantification
    call perform_kallisto_quantification {
      input:
        output_directory = output_directory_stripped,
        transcriptome_index = build_kallisto_index.reference_index,
        file_1 = unzip_file.file_1,
        file_2 = unzip_file.file_2,
        sample_name = sample_name
    }
  }

  # transform transcript counts to gene counts (to be used in edgeR) and also run differential expression
  # via Sleuth and edgeR
  call counts_and_differential_expression_output {
    input:
      DE_method = DE_method,
      output_directory = output_directory_stripped,
      transcript_counts_tar = perform_kallisto_quantification.transcript_counts_tar,
      organism = organism,
      samples_described_file = samples_described_file,
      logCPM_threshold = logCPM_threshold,
      percentage_aligned_threshold = percentage_aligned_threshold,
      ensemble_mirror = ensemble_mirror
  }

  # perform multiqc
  call perform_multiqc {
    input:
      output_directory = output_directory_stripped,
      transcript_counts_tar = perform_kallisto_quantification.transcript_counts_tar
  }

  output {
    String output_dir = output_directory_stripped
  }
}

task run_bcl_convert_and_define_read_pairs {
  input {
    File? bcl_sample_sheet
    String input_bcl_directory
    String output_directory
    Boolean no_lane_splitting
    Boolean strict_mode
    String bcl_convert_version
    String input_data_type
    String read_pairs_file

    String zones
    Int num_cpu
    String memory
    Int disk_space
    Int preemptible
  }

  String run_id = basename(input_bcl_directory)

  command <<<

    set -e
    export TMPDIR=/tmp

    if [[ '~{input_data_type}' == 'bcl' ]]; then

      gsutil -q -m cp -r ~{input_bcl_directory} .

      bcl-convert --bcl-input-directory ~{run_id} \
                  --output-directory fastq \
                  --sample-sheet  ~{bcl_sample_sheet} \
                  ~{true="--no-lane-splitting true" false=""  no_lane_splitting} \
                  ~{true="--strict-mode true" false=""  strict_mode}

      gsutil -q -m cp -r fastq/ ~{output_directory}/

      # create read_pairs file and copy to google bucket
      fastq_files=$(find ./fastq -name "*.fastq.gz" ! -name "Undetermined*" -type f -exec basename {} \; | sort) # get list of non-'Undetermined' fastq files and ensure they're paired together via sorting
      echo $fastq_files
      echo $fastq_files > read_pairs.tsv # write fastq names to file
      sed "s/ /\n/g" read_pairs.tsv  > read_pairs2.tsv # replace every whitespace with a newline character
      out_dir=$(echo '~{output_directory}' | sed 's;/;\\/;g' ) 
      sed "s/^/$out_dir\/fastq\//g" read_pairs2.tsv > read_pairs3.tsv  # add path to google bucket at beginning of each file name
      sed "N;s/\n/\t/" read_pairs3.tsv > read_pairs.tsv  # for every line corresponding to 'Read 1', replace newline with tab
      gsutil -m cp read_pairs.tsv ~{read_pairs_file} # write read_pairs.tsv file we generated to the specified read_pairs_file location

    else
      # localize read_pairs file
      gsutil -m cp ~{read_pairs_file} ./read_pairs.tsv
    fi

  >>>

  output {
    File read_pairs_tsv = "read_pairs.tsv"
  }
  
  runtime {
    preemptible: preemptible
    bootDiskSizeGb: 12
    disks: "local-disk ${disk_space} HDD"
    docker: "gcr.io/microbiome-xavier/bcl_convert:${bcl_convert_version}"
    cpu: num_cpu
    zones: zones
    memory: memory
  }

}

task build_kallisto_index {
  input {
    File reference_transcriptome
    Int cycle_number
  }

  String reference_transcriptome_basename = sub(basename(reference_transcriptome), ".fa.gz", "")

  command <<<
    set -e

    # list reference indices in google bucket
    indices=$(gsutil ls gs://rnaseq_reference_data/kallisto_indices/)

    # check if reference index already exists in google cloud bucket
    index_exists=0
    for index in $indices; do
      if [[ $index == *"~{reference_transcriptome_basename}_cycle~{cycle_number}.idx" ]]; then
          index_exists=1
          break
      fi
    done

    # if reference index already exists, then copy from google bucket -- otherwise, generate the reference index and copy to google bucket
    if [! -z "$indices" ] && [[ $index_exists == 1 ]]; then
      echo "reference index already exists"
      # localize reference index
      gsutil -m cp gs://rnaseq_reference_data/kallisto_indices/~{reference_transcriptome_basename}_cycle~{cycle_number}.idx ./
    else
      echo "generating reference index"
      # generate reference index
      kallisto index --index=~{reference_transcriptome_basename}_cycle~{cycle_number}.idx \
                      --kmer-size=~{cycle_number} \
                      ~{reference_transcriptome}

      # copy reference index to google bucket
      gsutil -m cp ~{reference_transcriptome_basename}_cycle~{cycle_number}.idx gs://rnaseq_reference_data/kallisto_indices/
    fi
  >>>

  output {
    File reference_index = "${reference_transcriptome_basename}_cycle${cycle_number}.idx"
  }

  runtime {
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    zones: "us-east1-b us-east1-c us-east1-d"
    cpu: 4
    memory: "16GB"
    preemptible: 2
    disks: "local-disk 80 HDD"
  }
}

task unzip_file {
  input {
    File zipped_file_1
    File zipped_file_2
  }

  String unzipped_basename_1 = sub(basename(zipped_file_1), ".gz", "")
  String unzipped_basename_2 = sub(basename(zipped_file_2), ".gz", "")

  command <<<
    set -e

    if [[ ~{zipped_file_1} == *.gz ]]; then
      # unzip files if they are zipped
      unpigz -c ~{zipped_file_1} > ~{unzipped_basename_1}
      unpigz -c ~{zipped_file_2} > ~{unzipped_basename_2}
    fi
  >>>

  output {
    File file_1 = unzipped_basename_1
    File file_2 = unzipped_basename_2
  }

  runtime {
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    zones: "us-east1-b us-east1-c us-east1-d"
    cpu: 4
    memory: "8GB"
    preemptible: 2
    disks: "local-disk 80 HDD"
  }
}

task perform_fastqc {
  input {
    String output_directory
    File file_1
    File file_2
    String sample_name
  }

  command <<<
    set -e

    # perform fastqc
    mkdir -p ~{sample_name}
    fastqc -f fastq -o ~{sample_name} ~{file_1} ~{file_2}

    # tar fastqc output directory
    tar zcfv fastqc_~{sample_name}.tar.gz ~{sample_name}

    # copy fastqc output to google bucket
    gsutil -q -m cp -r ~{sample_name} ~{output_directory}/fastqc/
  >>>

  output {
    File fastqc_output = "fastqc_${sample_name}.tar.gz"
  }

  runtime {
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    zones: "us-east1-b us-east1-c us-east1-d"
    cpu: 4
    memory: "4GB"
    preemptible: 2
    disks: "local-disk 80 HDD"
  }
}

task perform_kallisto_quantification {
  input {
    String output_directory
    File transcriptome_index
    File file_1
    File file_2
    String sample_name
  }

  command <<<
    set -e

    # perform quantification
    mkdir -p transcript_counts_~{sample_name}
    kallisto quant --index=~{transcriptome_index} \
                    --output-dir=transcript_counts_~{sample_name} \
                    --bootstrap-samples=100 \
                    --threads=11 \
                    ~{file_1} ~{file_2} > transcript_counts_~{sample_name}/stdout_qc.txt 2>&1

    # tar output directory
    tar zcfv transcript_counts_~{sample_name}.tar.gz transcript_counts_~{sample_name}

    # copy output to google bucket
    gsutil -q -m cp transcript_counts_~{sample_name}/* ~{output_directory}/transcript_counts/~{sample_name}/
  >>>

  output {
    File transcript_counts_tar = "transcript_counts_${sample_name}.tar.gz"
  }

  runtime {
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    zones: "us-east1-b us-east1-c us-east1-d"
    cpu: 12
    memory: "8GB"
    preemptible: 2
    disks: "local-disk 80 HDD"
  }
}

task counts_and_differential_expression_output {
  input {
    String DE_method
    String output_directory
    Array[File] transcript_counts_tar
    String organism
    String samples_described_file
    Float logCPM_threshold
    Float percentage_aligned_threshold
    String ensemble_mirror
  }

  command <<<
    set -e -o pipefail

    # NOTE: not yet sure how to use output from perform_kallisto_quantification directly here,
    #       doing so would probably be a better solution than copying from google bucket

    # copy transcript counts from google bucket
    # TODO: should probably change this to using output from last step
    gsutil -m cp -r ~{output_directory}/transcript_counts ./

    if [ ! -z "~{samples_described_file}" ]; then
      # copy samples described files
      gsutil -m cp ~{samples_described_file} ./samples_described.tsv
      echo "all_vs_all" > ./samples_compared.tsv
    else
      touch ./samples_described.tsv samples_compared.tsv
    fi

    # run kallisto_output
    mkdir -p /cromwell_root/kallisto_output
    Rscript /scripts/get_kallisto_counts.R ./transcript_counts ~{organism} samples_described.tsv samples_compared.tsv /cromwell_root/kallisto_output ~{ensemble_mirror}
    # tar output
    tar zcfv kallisto_output.tar.gz /cromwell_root/kallisto_output
    # copy output to google bucket
    gsutil -m cp -r /cromwell_root/kallisto_output ~{output_directory}/

    # run DE analysis .
    mkdir -p /cromwell_root/DE_output
    if [ -s samples_described.tsv ] && [ -s samples_compared.tsv ]; then
      Rscript /scripts/run_DE_analysis.R ~{DE_method} /cromwell_root/kallisto_output/kallisto_gene_counts_rounded.csv ./transcript_counts  ~{organism} samples_described.tsv samples_compared.tsv ~{logCPM_threshold} ~{percentage_aligned_threshold} /cromwell_root/DE_output
    fi

    # tar output
    tar zcfv DE_output.tar.gz /cromwell_root/DE_output
    # copy output to google bucket
    if [ -s samples_described.tsv ] && [ -s samples_compared.tsv ]; then
      gsutil -m cp -r /cromwell_root/DE_output ~{output_directory}/
    fi
  >>>

  output {
    File kallisto_counts = "kallisto_output.tar.gz"
    File DE_output = "DE_output.tar.gz"
  }

  runtime {
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    zones: "us-east1-b us-east1-c us-east1-d"
    cpu: 4
    memory: "100GB"
    preemptible: 2
    disks: "local-disk 100 HDD"
  }
}

task perform_multiqc {
  input {
    String output_directory
    Array[File] transcript_counts_tar
  }

  command <<<
    set -e

    mkdir -p ./multiqc_output

    # copy fastqc info
    gsutil -q -m cp -r ~{output_directory}/fastqc ./multiqc_output

    # copy kallisto info
    gsutil -q -m cp -r ~{output_directory}/transcript_counts ./multiqc_output

    # run multiqc
    cd ./multiqc_output
    multiqc ./
    cd ..

    # push multiqc report to google bucket
    gsutil -q -m cp ./multiqc_output/multiqc_report.html ~{output_directory}/
  >>>

  output {
    File multiqc_output = "multiqc_output/multiqc_report.html"
  }

  runtime {
    docker: "gcr.io/genomics-xavier/kco_bulk_rna_seq"
    zones: "us-east1-b us-east1-c us-east1-d"
    cpu: 4
    memory: "8GB"
    preemptible: 2
    disks: "local-disk 80 HDD"
  }
}
