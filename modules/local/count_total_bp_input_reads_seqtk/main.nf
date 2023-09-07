process COUNT_TOTAL_BP_INPUT_READS_SEQTK {

    publishDir   "${params.process_log_dir}",
        mode:    "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs:  { filename -> "${meta.id}.${task.process}${filename}" }

    tag { "${meta.id}" }

    container "gregorysprenger/seqtk@sha256:756bff7222c384d358cb22ecbbae443e112b296503cb0e1a6baf9cf80545ae20"

    input:
    tuple val(meta), path(reads), path(qc_input_filecheck)

    output:
    path ".command.out"
    path ".command.err"
    path "versions.yml"                             , emit: versions
    tuple val(meta), path("${meta.id}.total_bp.txt"), emit: total_bp

    shell:
    total_bp = 0
    '''
    source bash_functions.sh

    # Calculate total bp input for R1 FastQ file
    seqtk fqchk \
      !{reads[0]} \
      1> seqtk-fqchk.!{meta.id}.stdout.log \
      2> seqtk-fqchk.!{meta.id}.stderr.log

    # Extract just the total bp count of the R1 input FastQ file, then
    #  double it to estimate total R1 and R2 input
    if [ -s seqtk-fqchk.!{meta.id}.stdout.log ] ; then
      R1_total_bp=$(grep '^ALL' seqtk-fqchk.!{meta.id}.stdout.log | awk '{print $2}')
      if ! [[ $R1_total_bp =~ ^[0-9]+$ ]]; then
        msg "ERROR: R1 total bp size not counted with seqtk fqchk" >&2
        exit 1
      else
        # Double the R1 and skip R2 count; imperfect but faster
        #  and close enough for an estimation
        R1R2_total_bp=$(( 2 * ${R1_total_bp} ))
      fi
    else
      msg "ERROR: nucleotide count output logfile by seqtk fqchk is empty" >&2
      exit 1
    fi

    echo -n "${R1R2_total_bp}" > !{meta.id}.total_bp.txt

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        seqtk: $(seqtk 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
