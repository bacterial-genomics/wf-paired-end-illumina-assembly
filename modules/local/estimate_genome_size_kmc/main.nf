process ESTIMATE_GENOME_SIZE_KMC {

    tag { "${meta.id}" }
    container "gregorysprenger/kmc@sha256:27603041f8c8818aa71a1d0386df17eddca59dbd6441b7e84b78b8a09dc137df"
    label "process_medium"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.genome_size.txt"), emit: genome_size
    path(".command.{out,err}")
    path("versions.yml")                               , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Calculate unique kmers in the R1 read file as a proxy for genome size
    msg "INFO: Using kmc to estimate genome size"
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
        msg "ERROR: Genome size not estimated with kmc" >&2
        exit 1
      fi
    else
      msg "ERROR: kmer count output logfile by kmc is empty" >&2
      exit 1
    fi
    rmdir --ignore-fail-on-non-empty kmc-tmp-dir.!{meta.id}
    rm -rf kmc-binary-output-prefix.!{meta.id}*

    # Report the estimated genome size
    msg "INFO: Estimated genome size of !{meta.id}: ${genome_size}"

    echo -n "${genome_size}" > "!{meta.id}.genome_size.txt"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        kmc: $(kmc --help | grep "ver\\." | sed 's/.*ver\\. //1' | awk '{print $1}')
    END_VERSIONS
    '''
}
