# CheckM2

## Introduction

This process uses [CheckM2](https://github.com/chklovski/CheckM2) published in [2023](https://pubmed.ncbi.nlm.nih.gov/37500759/) to assess the full assembly file with a machine learning classfication algorithm for completeness and contamination. Assessment relies on the pre-computed diamond database, where protein fragments are used to train for contamination estimates with gradient boosting (GB) as well as for completeness estimates using neural networks (NN) and GB.

> "CheckM2 was far more accurate for medium, low-quality and highly contaminated genomes than[_sic_] both other [checkm and busco] tools"

## How CheckM2 works

From [CheckM2's documentation](https://github.com/chklovski/CheckM2):
> CheckM2 uses two distinct machine learning models to predict genome completeness.
>
> - The 'general' gradient boost model is able to generalize well and is intended to be used on organisms not well represented in GenBank or RefSeq (roughly, when an organism is novel at the level of order, class or phylum).
> - The 'specific' neural network model is more accurate when predicting completeness of organisms more closely related to the reference training set (roughly, when an organism belongs to a known species, genus or family).
> CheckM2 uses a cosine similarity calculation to automatically determine the appropriate completeness model for each input genome, but you can also force the use of a particular completeness model, or get the prediction outputs for both.
> There is only one contamination model (based on gradient boost) which is applied regardless of taxonomic novelty and works well across all cases.
