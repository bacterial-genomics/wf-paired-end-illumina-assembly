process CALCULATE_METRICS_FASTQ_SEQKIT {

    tag { "${meta.id}" }
    container "staphb/seqkit@sha256:8eb09a52ae932f7c25cfbb8db0df7110567087a187c7e90d46f499962d1c82c9"

    input:
    tuple val(meta), path(reads)
    val(input_fastq_type)

    output:
    path("${meta.id}.${input_fastq_type}.metrics_summary.tsv"), emit: output
    path(".command.{out,err}")
    path("versions.yml")                                      , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Calculating !{meta.id} statistics of !{input_fastq_type} FastQ input with SeqKit..."

    fastq_files=( !{reads} )
    num_reads="${#fastq_files[@]}"

    msg "INFO: Found ${num_reads} of !{input_fastq_type} FastQ files: ${fastq_files[@]}"
    msg "INFO: Found !{input_fastq_type} FastQ files: !{reads}"

    # Calculate stats on 1 or more FastQ input files
    # SeqKit prints each file stats separately, line-by-line, no total.
    # Unlike seqtk, SeqKit autohandles input files with or without compression
    #   and reports correct numbers regardless (thank you! @shenwei356)
    seqkit \
      stats \
      --tabular \
      --threads "!{task.cpus}" \
      --out-file "!{meta.id}.!{input_fastq_type}.seqkit_stats.tsv" \
      !{reads}

    msg "INFO: Calculated statistics with SeqKit for !{input_fastq_type} FastQ files: !{reads}"
    msg "INFO: Calculated statistics with SeqKit for ${num_reads} !{input_fastq_type} FastQ files: ${fastq_files[@]}"

    awk -v sample_id="!{meta.id}" '
    BEGIN {
        # Set the header contents (renamed)
        OFS = "\t"
        header_translation["num_seqs"] = "Total_Length_[bp]"
        header_translation["sum_len"]  = "Total_Sequences_[#]"
        header_translation["min_len"]  = "Minimum_Sequence_Length_[bp]"
        header_translation["avg_len"]  = "Mean_Sequence_Length_[bp]"
        header_translation["max_len"]  = "Maximum_Sequence_Length_[bp]"
    }

    NR == 1 {
        # Change header row item from "file" into "Sample_name"
        $1 = "Sample_name"  

        # Rename specific header names
        for (i = 1; i <= NF; i++) {
            if ($i in header_translation) {
                $i = header_translation[$i]
            }
        }

        # Print the modified header, excluding columns 2 and 3 ("format" and "type")
        print $1, $4, $5, $6, $7, $8
    }

    NR > 1 {
        # Process data rows, excluding columns 2 and 3 ("format" and "type")
        # Aggregate results if there are multiple rows (from >1 FastQ input)
        num_seqs += $4
        sum_len  += $5

        if (NR == 2) {
            # First row: initialize min/max/avg values
            min_len = $6
            max_len = $8
            avg_len = $7
        } else {
            # Update min and max
            if ($6 < min_len) min_len = $6
            if ($8 > max_len) max_len = $8
        }
    }

    END {
        # Always print the aggregated results (or single line if only one row)
        print sample_id, sum_len, num_seqs, min_len, avg_len, max_len
    }' \
    "!{meta.id}.!{input_fastq_type}.seqkit_stats.tsv" \
    > "!{meta.id}.!{input_fastq_type}.metrics_summary.tsv"

    filepath="$(readlink -f !{meta.id}.!{input_fastq_type}.metrics_summary.tsv)"
    msg "INFO: Summarized SeqKit statistics of !{meta.id} !{input_fastq_type} for ${filepath}"

    # NOTE: DONE!
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
        seqkit: $(seqkit 2>&1 | grep "^Version: " | sed 's/^Version: //1')
    END_VERSIONS
    '''
}
