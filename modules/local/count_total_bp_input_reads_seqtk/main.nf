process COUNT_TOTAL_BP_INPUT_READS_SEQTK {

    tag { "${meta.id}" }
    container "staphb/seqtk@sha256:82797114adb664ba939b4f6dfcb822483a4af827def4288e5207be559055f2cc"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.input_total_bp.txt"), emit: input_total_bp
    path(".command.{out,err}")
    path("versions.yml")                                  , emit: versions

    shell:
    total_bp = 0
    '''
    source bash_functions.sh

    # Calculate total bp input for R1 FastQ file
    # NOTE: specified quality to trim param depends on user-supplied read
    #       trimmer selected (e.g., Fastp or Trimmomatic) and each have
    #       different variables for minimum quality score to use:
    #           - Fastp       = params.fastp_window_mean_quality
    #           - Trimmomatic = params.trimmomatic_required_quality
    #       but we'll just use Phred30 cutoff to be conservative for the
    #       downsampling step here (`-q 30`).

    msg "INFO: Calculating !{meta.id} basepairs above Phred 30 with Seqtk for subsampling calculations..."

    seqtk fqchk \
      -q 30 \
      !{reads[0]} !{reads[1]} \
      > seqtk-fqchk.!{meta.id}.stdout.log \

    msg "INFO: Calculated !{meta.id} basepairs above Phred 30 with Seqtk"

    # Extract just the total bp count of the R1 input FastQ file, then
    #  double it to estimate total R1 and R2 input
    if [ -s seqtk-fqchk.!{meta.id}.stdout.log ] ; then
      total_bp=$(grep '^ALL' seqtk-fqchk.!{meta.id}.stdout.log | awk '{printf $2}')
      if ! [[ $total_bp =~ ^[0-9]+$ ]]; then
        msg "ERROR: total bp = $total_bp" >&2
        msg "ERROR: total bp size not counted with seqtk fqchk" >&2
        exit 1
      else
        echo -n "${total_bp}" > "!{meta.id}.input_total_bp.txt"
        msg "INFO: found ${total_bp}bp for !{meta.id}"
      fi
    else
      msg "ERROR: nucleotide count output logfile by seqtk fqchk is empty" >&2
      exit 1
    fi

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        seqtk: $(seqtk 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
