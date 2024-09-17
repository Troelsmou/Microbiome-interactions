# HAPk

This is the code for calculating the HAPk heuristic.

## Software Requirements

This script was tested with:  
R==4.4.1  
tidyverse==2.0.0

## Installation

Clone this repository.  
Install R and Rstudio.  
Install the required packages from cran.  
And then you can run the scripts.

## Usage

```bash
Rscript --args "20220330_lotus_snps.csv" 100 "Reproduced" 8
```

## Input files

`20220330_lotus_snps.csv` [Numeric SNP calls (on figshare)]

## Output files

`Reproduced_cluster_result.csv` [Csv file with 100 kb windowed clustering R-squared values for 2-8 clusters and HAPk value]

## Runtime
Typical runtime: 40 minutes
