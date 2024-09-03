#!/usr/bin/env bash

summarize_kraken() {
  # $1 = input tab-delimited kraken data file (e.g., "*.kraken2_output.tsv") from `kraken --report <OUTPUT TSV>`

  # Initialize an array to collect all output sections
  output=()

  # Report unclassified
  UNCL=($(grep $'\tU\t' "${1}" | head -n1 | awk '{print $1,$2,$6}'))
  if [[ ${#UNCL[@]} -eq 3 ]]; then
    tabline=$(echo "${UNCL[@]}" | sed -E 's/ +/%\t/1' | sed -E 's/ +/\t/1')
    output+=("$tabline")
  else
    output+=("0%\t0\tUnclassified")
  fi

  # At most top 3 genera
  while read -r l; do
    tabline=$(echo "${l}" | sed -E 's/ +/%\t/1' | sed -E 's/ +/\t/1')
    output+=("$tabline")
  done < <(grep $'\tG\t' "${1}" | head -n3 | awk -F $'\t' '{print $1,$2,$6}')

  # At most top 3 species
  while read -r l; do
    tabline=$(echo "${l}" | sed -E 's/ +/%\t/1' | sed -E 's/ +/\t/1')
    output+=("$tabline")
  done < <(grep $'\tS\t' "${1}" | head -n3 | awk -F $'\t' '{print $1,$2,$6}')

  # Join the array with tabs and print the final result
  (
    IFS=$'\t'
    echo "${output[*]}"
  )
}
