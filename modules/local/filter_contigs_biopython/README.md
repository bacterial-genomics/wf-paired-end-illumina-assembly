# Contig filtering using biopython script

## Introduction

Genome assemblers often save all contiguous sequences that form (from 1 or more overlap). Since the FastQ data of a bacterial isolate are generally large enough that sequencing errors ("noise") do show up, small levels of unrelated taxa can get mixed in ("contamination"), and a variety of others issues occur, it is almost never a good bet to take direct assembler output and use it.

I have written a standalone biopython script [here](https://github.com/chrisgulvik/genomics_scripts/blob/main/filter.contigs.py) to evaluate each contig for compositional complexity (requiring at least 3 nucleotides to be present), excessive GC skew, a minimum 1 kbp length, and at least 5x coverage.

The reason for each contig discarded is logged, and a user can even store all discarded contigs in a separate FastA file to evaluate each contig record individually. A variety of calculated statistics are also determined, which enables rapid assessment on removed and retained contigs.

## Log information

- Input number and lengths of all contigs
- Discarded number and lengths of all contigs
- Output number and lengths of all contigs
- Which of the 4 popular assembler contig formats (IDBA, SKESA, Unicycler, Velvet/SPAdes) was automatically detected for defline data parsing
- The reason why each contig discarded (e.g., "too short")
- Coverage quartiles (25, 50, 75) of all contigs that passed all filters
- The minimum, average, and maximum depths observed in all contigs that passed all filters
