# Host removal using Hostile

## About this module

This process uses the [Hostile](https://github.com/bede/hostile) package to align FastQ reads to a reference genome and stores all unmapped reads. Its default use is to remove human reads, with the primary output being FastQ files.

## Reference Genome

- One terrific feature of this package is its bundled assembly files of human genomes, and the default uses the latest telomere-to-telomere (T2T) human assembly that Adam Phillippy and colleagues accomplished in [2022](https://pubmed.ncbi.nlm.nih.gov/35357919/) combined with human
  leukocyte antigen (HLA) sequences from [here](https://www.ebi.ac.uk/ipd/imgt/hla/release/v351/).

- If a user has a different background genome to remove, for example, perhaps the bacterial pathogen was collected from a brown common sewer rat from New York City (_Rattus norvegicus_) or a Holstein dairy cow in the USA (_Bos taurus_), a non-default reference genome can be supplied. For Illumina read removal with a non-default reference genome, a user must specify a path prefix of all 6 of the bowtie2 index files-- not a FastA file.

- To avoid removing reads that match to both the target pathogen and the host and maximize retention of the target pathogen, `hostile` supplies the utilities to compare both at the same time and form a custom reference file [here](https://github.com/bede/hostile#masking-reference-genomes).

A clever approach in Bede Constantinides preprint paper was to gather [985](https://github.com/bede/hostile/blob/main/paper/supplementary-table-2.tsv) reference grade bacterial genomes from [the FDA ARGOS database](https://www.ncbi.nlm.nih.gov/bioproject/231221) into a large FastA for masking the combined T2T-CHM and IPD-IMGT/HLA human genome reference, which is available [here](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985.tar).

## Software Citation

A peer-reviewed manuscript is likely in progress, but for now the [biorxiv paper](https://www.biorxiv.org/content/10.1101/2023.07.04.547735v1) describes the Hostile package.
