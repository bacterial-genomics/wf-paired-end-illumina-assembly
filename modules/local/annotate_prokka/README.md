# Genome annotation using Prokka

## Introduction

This process uses [prokka](https://github.com/tseemann/prokka) published in [2014](https://pubmed.ncbi.nlm.nih.gov/24642063/) to annotate a FastA assembly file. It identifies coding sequences (CDS), ribosomal RNAs (rRNA), transfer RNAs (tRNA), and other features. The primary output is a GenBank formatted file.

## Alignment/similarity scoring

To boost annotation confidence, we default to an 1e-08 as a similarity e-value cut-off, but this parameter can be modified if a user wishes to relax the homology assignments which could be useful for novel taxa or highly distant homolog matching.

## Alternatives

Prokka is meant for rapid annotation, which performs very well in most cases. However, the NCBI's annotation pipeline "Prokaryotic Genome Annotation Pipeline" [(PGAP)](https://github.com/ncbi/pgap) is recommended for the most comprehensive annotations, especially if you believe the functional name is to a distant homolog or if the gene of interest appears to be missing in Bakta's GenBank output file.

PGAP was first developed in 2001 by NCBI and [GATech](https://www.gatech.edu/) staff and published in [2016](https://pubmed.ncbi.nlm.nih.gov/27342282/). The most recent updates have been reported in [2021](https://pubmed.ncbi.nlm.nih.gov/33270901/).
