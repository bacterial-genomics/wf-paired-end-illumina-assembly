process CALCULATE_COVERAGE_UNIX {

    tag { "${meta.id}-${meta.assembler}" }
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(summary_assembly_metrics), path(summary_cleanedreads), path(summary_aligned_stats)

    output:
    path("${meta.id}-${meta.assembler}.GenomeCoverage.tsv"), emit: summary
    path(".command.{out,err}")
    path("versions.yml")                                   , emit: versions

    shell:
    '''
    source bash_functions.sh

    # Use assembly report sample names
    echo -n '' > "!{meta.id}-!{meta.assembler}.GenomeCoverage.tsv"
    i=0
    while IFS=$'\t' read -r -a ln; do
      msg "INFO: looking for depth of coverage for ${ln[0]}"

      # Skip unusually named assembly filenames that shouldn't be in this workflow
      if grep -q -e "unicyc_" -e ".uncorrected" <<< "${ln[0]}"; then
        msg "INFO: skipping ${ln[0]} due to containing unicyc_ or .uncorrected"
        continue
      fi

      # If the mapped/aligned stats file exists, use those depth values
      mapped_coverage=$(grep ${ln[0]} !{summary_aligned_stats} 2> /dev/null | awk '{print $9"x_+/-_" $10"x"}')

      if [[ "${mapped_coverage}" =~ ^[0-9]+[.][0-9]{1}x_[+][/][-]_[0-9]+[.][0-9]{1}x$ ]]; then
        msg "INFO: Cleaned read sequences mapped to assembly Coverage data found for !{meta.id}-!{meta.assembler}: ${mapped_coverage}"
        echo -e "${ln[0]}\t${mapped_coverage}" >> "!{meta.id}-!{meta.assembler}.GenomeCoverage.tsv"
        ((i=i+1))

      # If mapped/aligned stats are missing, revert to just extracting the
      #   cleaned reads basepair count and the assembly cumulative length to
      #   estimate coverage.
      # NOTE: This shouldn't happen, but it's a backup and easily discernible
      #       from the mapped-based method, because no stdev ("_+/-_") is
      #       reported with this method.
      else
        msg "INFO: Read alignment data absent for ${ln[0]}, so cleaned bases"
        msg "      given to the assembler were used to calculate coverage"
        basepairs=$(grep ${ln[0]} !{summary_cleanedreads} | cut -f 2)
        assembly_length=$(grep ${ln[0]} !{summary_assembly_metrics} | cut -f 4)
        if [[ -z "$basepairs" || -z "$assembly_length" || ! "$basepairs" =~ ^[0-9]+$ || ! "$assembly_length" =~ ^[0-9]+$ || "$basepairs" -le 0 || "$assembly_length" -le 0 ]]; then
          msg "ERROR: skipping ${ln[0]}: $basepairs bp or $assembly_length bp are unset, empty, not integers, or not greater than zero" >&2
          continue
        fi
        cleaned_basepairs_per_assembly_length=$(awk -v a="${basepairs}" -v b="${assembly_length}" '{print a / b}')
        cleaned_basepairs_per_assembly_length=$(printf "%.1fx" "$cleaned_basepairs_per_assembly_length")

        if [[ "${cleaned_basepairs_per_assembly_length}" =~ ^[0-9]+[.][0-9]{1}x$ ]]; then
          msg "INFO: Cleaned read basepairs per assembly site Coverage for !{meta.id}-!{meta.assembler}: ${cleaned_basepairs_per_assembly_length}"
          echo -e "${ln[0]}\t${cleaned_basepairs_per_assembly_length}" >> "!{meta.id}-!{meta.assembler}.GenomeCoverage.tsv"
          ((i=i+1))
        fi
      fi
    done < <(grep -v -e 'Total_length' -e 'Total length' !{summary_assembly_metrics})

    msg "INFO: stored coverage values for ${i} sample(s)"

    # Add header row to data output
    sed -i '1i Sample_name\tCoverage_(mean[x]_+/-_stdev[x])' "!{meta.id}-!{meta.assembler}.GenomeCoverage.tsv"

    # Get process version information
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
