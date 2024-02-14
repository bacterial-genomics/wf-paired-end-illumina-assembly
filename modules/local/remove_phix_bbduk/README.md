# PhiX read removal using BBDuk

## Introduction

The bacteriophage PhiX174 is a commonly added spike-in sequence for Illumina sequencing. When added, it serves as a positive control for the whole sequencing run. Many Illumina instruments align reads against this well-known short sequence (thanks to [Fred Sanger](https://pubmed.ncbi.nlm.nih.gov/731693/))to identify SNPs as a proxy for how reliable the unknown samples sequences are.

At one point, reported [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4511556/), the forgetful removal of this resulted in >1,000 genomes contaminated with it on NCBI and 10% of the genomes published in the literature. This artificially added sequence, when retained, will form overlaps of reads during assembly and join biological sample DNA that are not contiguous but would otherwise seem to be from the raw data. Therefore PhiX fragments must be removed to form a higher quality assembly.

## PhiX removal

Sequencing errors can (and do) occur occasionally, so I allow for a 1 SNP to occur in a 31 bp aligned stretch to PhiX.

- BBDuk here scans each FastQ read sequence for 31-mers that match to PhiX. When an exact or 1 SNP match occurs, the aligned portion is removed from the sequence, but the remaining unaligned component of the sequence is retained.
- I use this FastA file for the PhiX reference [here](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1/), which is 5,386 bp long and is the only assembly on NCBI (of ~150) that is the International Committee on Taxonomy of Viruses species exemplar.

### Log information

- Total input read count
- Total input base pair count
- Total base pairs identified as PhiX and discarded
