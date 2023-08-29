process ESTIMATE_GENOME_SIZE_KMC {

    publishDir   "${params.process_log_dir}",
        mode:    "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs:  { filename -> "${meta.id}.${task.process}${filename}" }

    tag { "${meta.id}" }

    container "nanozoo/kmc@sha256:e6375ae53d453cd9d3b3c0ce0890cf2c3fd259b8d459e4c081b2737d4b34979f"

    input:
    tuple val(meta), path(reads), path(qc_input_filecheck)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                                , emit: versions
    tuple val(meta), path("${meta.id}.genome_size.txt"), emit: genome_size

    shell:
    '''
    source bash_functions.sh

    # Calculate unique kmers in the R1 read file as a proxy for genome size
    mkdir -p kmc-tmp-dir.!{meta.id}
    kmc \
      -t!{task.cpus} \
      "!{reads[0]}" \
      kmc-binary-output-prefix.!{meta.id} \
      kmc-tmp-dir.!{meta.id} \
      1> kmc.!{meta.id}.stdout.log \
      2> kmc.!{meta.id}.stderr.log

    # Extract just the integer of unique kmers identified in the R1 FastQ
    if [ -s kmc.!{meta.id}.stdout.log ] ; then
      genome_size=$(grep 'unique counted k' kmc.!{meta.id}.stdout.log \
        | cut -d ':' -f 2 \
        | sed 's/[[:space:]]//g')

      if ! [[ ${genome_size} =~ ^[0-9]+$ ]]; then
        msg "ERROR: genome size not estimated with kmc" >&2
        exit 1
      fi
    else
      msg "ERROR: kmer count output logfile by kmc is empty" >&2
      exit 1
    fi
    rmdir --ignore-fail-on-non-empty kmc-tmp-dir.!{meta.id}
    rm -rf kmc-binary-output-prefix.!{meta.id}*

    # Report the estimated genome size
    msg "INFO: Genome size was estimated to be ${genome_size} bp with kmc"

    echo -n "${genome_size}" > !{meta.id}.genome_size.txt

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        kmc: $(kmc --help | grep "ver\\." | sed 's/.*ver\\. //1' | awk '{print $1}')
    END_VERSIONS
    '''
}
