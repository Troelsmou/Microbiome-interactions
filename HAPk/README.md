# Mantel-test

This is the code for doing the mantel test for the paper.

## Software Requirements

The needed packages are listed in requirements.txt

## Installation

Clone this repository.  
Install R and Rstudio.  
Install the packages in requirements.txt from cran.  
And then you can run the scripts.

## Usage

The three scripts are supposed to be run interactively in an R session.

## Input files

`Lotus_accessions_location_2.csv` [Lotus accessions and their sampled longitude and latitude]  
`lotusPop.csv` [Lotus accessions and their population membership]  
`microbiomeSNPsExtended.csv` [Filtered microbiome gwas results]  
`SNPs_to_generate_random_sets.txt` [List of SNPs fulfilling GWAS filtering criteria (close to or in a gene)]  
`20220330_lotus_snps.csv` [Numeric SNP calls (on figshare)]

## Output files

`fig_2j.pdf` [Histogram of mantel test null distribution with GWAS mantel score plotted as vertical line]

## File relationships

![PlantUML diagram](http://www.plantuml.com/plantuml/proxy?cache=no&src=https://raw.githubusercontent.com/Troelsmou/Microbiome-interactions/main/Mantel-test/diagram.puml)

## Runtime
Typical runtime: 10 minutes
