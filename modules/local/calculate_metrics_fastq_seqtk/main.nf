process CALCULATE_METRICS_FASTQ_SEQTK {

    tag { "${meta.id}" }
    container "staphb/seqtk@sha256:82797114adb664ba939b4f6dfcb822483a4af827def4288e5207be559055f2cc"

    input:
    tuple val(meta), path(reads)
    val(input_fastq_type)

    output:
    path("${meta.id}.${input_fastq_type}.seqtk.metrics_summary.tsv"), emit: output
    path(".command.{out,err}")
    path("versions.yml")                                            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Calculating !{meta.id} statistics of !{input_fastq_type} FastQ input with Seqtk..."

    fastq_files=( !{reads} )
    num_reads="${#fastq_files[@]}"

    msg "INFO: Found ${num_reads} of !{input_fastq_type} FastQ files: ${fastq_files[@]}"
    msg "INFO: Found !{input_fastq_type} FastQ files: !{reads}"

    # Calculate stats on 1 or more FastQ input files
    # Seqtk prints only totals regardless of 1, 2, or 3 input FastQ files, and
    #   gives 0-code exit status for files with GZ input but wrong numbers
    #   unless you pipe them in from stdin.
    if [[ !{reads[0]} == *.gz ]]; then
        zcat !{reads} | seqtk size > "!{meta.id}.!{input_fastq_type}.Seqtk_stats.tsv"
    else
        seqtk size !{reads} > "!{meta.id}.!{input_fastq_type}.Seqtk_stats.tsv"
    fi

    msg "INFO: Calculated statistics with Seqtk for !{input_fastq_type} FastQ files: !{reads}"

    awk -v sample_name="!{meta.id}" '
    BEGIN {
        # Add a header row of the FastQ sequence count data
        print "Sample_name\tTotal_Length_[bp]\tTotal_Sequences_[#]"
    }

    {
        # Print the data row
        print sample_name, $2, $1
    }' OFS="\t" \
    "!{meta.id}.!{input_fastq_type}.Seqtk_stats.tsv" \
    > "!{meta.id}.!{input_fastq_type}.seqtk.metrics_summary.tsv"

    filepath="$(readlink -f !{meta.id}.!{input_fastq_type}.seqtk.metrics_summary.tsv)"
    msg "INFO: Summarized Seqtk statistics of !{meta.id} !{input_fastq_type} for ${filepath}"

    # NOTE: DONE! This section can be removed later but keeping for now
    #       as a reminder of what this is supposed to accomplish.
    # # TO-DO: move this unix-only component to separate QA_READS_BASEPAIR_COUNT_UNIX
    # # Count nucleotides per read set
    # echo -n '' > "!{meta.id}-!{meta.assembler}.Clean_Reads-Bases.tsv"
    # for (( i=0; i<3; i+=3 )); do
    #   R1=$(basename "!{meta.id}_R1.paired.fq.gz" _R1.paired.fq.gz)
    #   R2=$(basename "!{meta.id}_R2.paired.fq.gz" _R2.paired.fq.gz)
    #   single=$(basename "!{meta.id}_single.fq.gz" _single.fq.gz)

    #   # Verify each set of reads groups properly
    #   nr_uniq_str=$(echo -e "${R1}\\n${R2}\\n${single}" | sort -u | wc -l)
    #   if [ "${nr_uniq_str}" -ne 1 ]; then
    #     msg "ERROR: improperly grouped ${R1} ${R2} ${single}" >&2
    #     exit 1
    #   fi
    #   echo -ne "${R1}\t" >> "!{meta.id}-!{meta.assembler}.Clean_Reads-Bases.tsv"
    #   zcat "!{meta.id}_R1.paired.fq.gz" "!{meta.id}_R2.paired.fq.gz" "!{meta.id}_single.fq.gz" | \
    #     awk 'BEGIN{SUM=0} {if(NR%4==2){SUM+=length($0)}} END{OFMT="%f"; print SUM}' \
    #       >> "!{meta.id}-!{meta.assembler}.Clean_Reads-Bases.tsv"

    #   sed -i '1i Sample_name\tCleaned_bases_(#)' "!{meta.id}-!{meta.assembler}.Clean_Reads-Bases.tsv"
    # done

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        seqtk: $(seqtk 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
