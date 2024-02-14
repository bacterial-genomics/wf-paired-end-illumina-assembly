# Adapter clipping and quality trimming using Trimmomatic

## Introduction

This process uses Trimmomatic to perform some of the read cleaning steps to remove and trim FastQ sequences.

## Adapter clipping

Trimmomatic supplies a multi-record FastA sequence of known Illumina adapter sequences [here](https://github.com/usadellab/Trimmomatic/tree/main/adapters) and some have concatenated all of those into a single adapter FastA file to identify and remove sequence reads containing any known adapters. The sequences they provide are just the adapter without barcodes. However, they do not supply sequences for some other kits that I have dealt with (e.g., Rubicon Genomics ThruPLEX, NEBNext).

I used NCBI's UniVec database [here](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) which contains Illumina adapter sequences and other unrelated sequences to form a more comprehensive multi-FastA adapter file, which includes the exact barcode names as well in their deflines for identifying which specific adapter(s) were removed. [This file](https://github.com/chrisgulvik/genomics_scripts/blob/master/examples/adapters_Nextera_NEB_TruSeq_NuGEN_ThruPLEX.fas.gz) is used for adapter removal.

The Illumina instruments should be detecting and removing perfect matches, but when there is a sequencing error or two, it can end up in the FastQ output. So, this Illumina adapter identification and removal process allows for up to 2 mismatches from each (roughly 50-70 bp lengths) sequence.

## Quality trimming

- A sliding window of 6 bp starts at the 5' sequence and requires a mean Phred score of 30 or higher to occur for it to be retained, and once the window moving to the 3' sequence hits a spot less than Phred 30, the remaining 3' sequence is discarded.
- If the first base pair of the trimmed read is below Phred 10, the sequence is removed.
- If the last base pair of the trimmed read is below Phred 10, the sequence is removed.

## Length filter

- After adapter clipping and quality trimming are performed, sequence reads might be so small that they are practically useless in an otherwise good set of sequences. The extra compute time normally will not improve the assembly, so trimmed reads that pass all the above filters less than 50 bp are discarded.

## Sister read pairing

- During the read cleaning processes Trimmomatic performs, sometimes 1 of the 2 reads in a pair were discarded whereas the other sister read passed all filters. This broken sister reads (or "singletons") are retained and passed along to the next process.

## Log information

Each of these trimming effects/outcomes are tallied and logged:

- Total discarded read count
- Number of forward reads that lack a high quality R2 sister read
- Number of reverse reads that lack a high quality R1 sister read
- Total number of broken read pairs saved as singletons

> [!NOTE]
> Trimmomatic does not list which specific adapters (by name) were removed but there is a [feature request](https://github.com/usadellab/Trimmomatic/issues/9) for this to be implemented within Trimmomatic.
