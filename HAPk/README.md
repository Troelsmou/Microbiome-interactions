# HAPk calculation

This is the code for calculating the HAPk heuristic.

## Software Requirements

The needed packages are listed in requirements.txt

## Installation

Clone this repository.  
Install R and Rstudio.  
Install the packages in requirements.txt from cran.  
And then you can run the scripts.

## Usage

```bash
Rscript --args "20220330_lotus_snps.csv" 100 "Reproduced" 8
```

## Input files

`20220330_lotus_snps.csv` [Numeric SNP calls (on figshare)]

## Output files

`Reproduced_cluster_result.csv` [Csv file with 100 kb windowed clustering R-squared values for 2-8 clusters and HAPk value]

## File relationships

![PlantUML diagram](http://www.plantuml.com/plantuml/proxy?cache=no&src=https://raw.githubusercontent.com/Troelsmou/Microbiome-interactions/main/Mantel-test/diagram.puml)

## Runtime
Typical runtime: 40 minutes
