# GWAS-filtering

This is the code for filtering the raw GWAS output to the curated results presented in the paper.

## Software Requirements

R==4.4.1  
tidyverse==2.0.0

## Installation

Clone this repository.  
Install R and Rstudio.  
Install the required packages.  
And then you can run the script.

## Usage

The script is supposed to be run interactively in an R session.

## Input files

`20210713_Lj_Gifu_v1.3_predictedGenes.gff3` [Lotus japonicus gff file containing gene locations (can be downloaded from https://lotus.au.dk/data/download)]  
`20240508_Lotus_rarefied_GWAS_Results_sig.csv.gz` [GWAS results on bacterial and non-bacterial traits, filtered to a p.value of 10^-4 or lower (on figshare)]

## Output files

`20240508_Lotus_rarefied_GWAS_Results_sig_filtered_bacteria.csv` [CSV file containing curated bacterial GWAS results]  
`20240508_Lotus_rarefied_GWAS_Results_sig_filtered_other.csv` [CSV file containing curated GWAS results from non-bacterial traits]


## Runtime
Typical runtime: 2 minutes
