# GWAS-filtering

This is the code for filtering the raw GWAS output to the curated results presented in the paper.

## Software Requirements

This script was tested with:  
R==4.4.1  
tidyverse==2.0.0
data.table==1.18.4

## Installation

Clone this repository.  
Install R and Rstudio.  
Install the required packages from cran.  
And then you can run the script.

## Usage

The script is supposed to be run interactively in an R session.

## Input files

`Lotus_GD.csv` [File containing all SNP calls for the 155 accessions (on zenodo)]
`Lotus_GM.csv` [File containing SNP chromosome and position for all SNPs (on zenodo)]
`20210713_Lj_Gifu_v1.3_predictedGenes.gff3` [Lotus japonicus gff file containing gene locations (can be downloaded from https://lotus.au.dk/data/download)]
`20240508_Lotus_rarefied_GWAS_Results_sig.csv` [GWAS results on bacterial and non-bacterial traits, filtered to a p.value of 10^-4 or lower (on zenodo)]
`20250408_Lotus_permutation_results.csv` [Permutation GWAS results on bacterial and non-bacterial traits (on zenodo)]
`20250714_Lotus_permutation_GWAS_Results_sig.csv` [GWAS results on permuted bacterial and non-bacterial traits, filtered to a p.value of 10^-4 or lower (on zenodo)]
`20250808_Lotus_nonsense_GWAS_permutation_results.csv` [Permutation GWAS results on permuted bacterial and non-bacterial traits (on zenodo)]
`Lotus_Johan_ave_rep_family_mds_rarefied.csv` [Phenotype file for both abiotic and biotic phenotypes across accessions (on zenodo)]

## Output files

`20250825_Permutation_filtered_nonsense_significants.csv` [CSV file containing curated permuted bacterial GWAS results]
`20250902_Permutation_filtered_nonsense_abiotic_significants.csv` [CSV file containing curated permuted abiotic GWAS results]
`20250514_Permutation_filtered_significants.csv` [CSV file containing curated bacterial GWAS results]
`20250514_Permutation_filtered_significants_non_bact.csv` [CSV file containing curated GWAS results from abiotic traits]


## Runtime
Typical runtime: 10 minutes
