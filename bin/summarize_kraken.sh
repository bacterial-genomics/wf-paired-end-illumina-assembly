#!/usr/bin/env bash

summarize_kraken()
{
	#BASE=$(basename $1 _kraken.tab)
	#echo -en "$BASE\t"

	# report unclassified
	UNCL=(`grep $'\tU\t' ${1} | head -n1 | awk '{print $1,$2,$6}'`)
	if [[ ${#UNCL[@]} -eq 3 ]]; then
		tabline=$(echo ${UNCL[@]} | sed -E 's/ +/%\t/1' | sed -E 's/ +/\t/1')
		echo -en "$tabline\t"
	else
		echo -en "0%\t0\tUnclassified\t"
	fi

	# At most top 3 genera
	while read -r l; do
		tabline=$(echo "${l}" | sed -E 's/ +/%\t/1' | sed -E 's/ +/\t/1')
		echo -en "$tabline\t"
	done < <(grep $'\tG\t' ${1} | head -n3 | awk -F $'\t' '{print $1,$2,$6}')

	# At most top 3 species
	while read -r l; do
		tabline=$(echo "${l}" | sed -E 's/ +/%\t/1' | sed -E 's/ +/\t/1')
		echo -en "$tabline\t"
	done < <(grep $'\tS\t' ${1} | head -n3 | awk -F $'\t' '{print $1,$2,$6}')
	echo ''
}
