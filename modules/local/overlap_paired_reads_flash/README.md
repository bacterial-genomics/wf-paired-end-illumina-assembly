# Merge overlapping reads using FLASH

## Introduction

During the library preparation of genomic DNA fragments, an important step involves isolating a size normally around 500 bp depending on the methods used. Generally, fragments too long have a lower clustering efficiency, quality (fidelity), and yield, which are not ideal [ref](https://pubmed.ncbi.nlm.nih.gov/30814542/).

When fragments are too short, during paired-end sequencing, the sister reads can end of having some of the same content due to sequencing overlap. Some assemblers and read mappers might be able to handle some of this, but ideally paired reads have a missing unknown sequence between the two sister reads.

When a pair has overlapping content, this module assembles (or overlaps) those pairs into longer singleton reads. Given the redundancy nature of bacterial genome content, the majority (80% or more) of the read pair must overlap for it to be collapsed into a singleton read.

## Log information

- read length (in base pairs) detected from first 100 raw input sequences in the forward (R1) FastQ
- Total pairs overlapped into singleton reads
- Total non-overlapping pairs
