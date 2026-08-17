## Prediction

This is the code for performing the genomic prediction gBLUP-style of traits used in the manuscript as a more robust way to ascertain the variance that is attributable to genetics in each trait. The script has an additional capability to predict a trait based on a set of SNPs from a different file, comparing it to random sets of SNPs equal in number to the set of SNPs from the file, but this capability is not used in the manuscript.

## Software Requirements

This script was developed and run with the following in mind:<br>
R==4.3.2<br>
tidyverse==2.0.0<br>
data.table==1.14.10<br>
pacman==0.5.1<br>
rrBLUP==4.6.3<br>
`
## Inputs

`Empty_sigfile.csv` [Header only CSV file, only contains data when the additional SNP-set capability is used]<br>
`20260615_Lotus_ave_rep_family_iso_mds_rarefied_newformat.csv` [Phenotype file of traits to predict]<br>
`Lotus_GD.csv` [File containing all SNP calls for the 155 accessions (on zenodo)]<br>

## Arguments

--pheno: The phenotype file<br>
--prediction_trait: Which trait to predict in the pheno file<br>
--GD: The GD-formatted genotype file<br>
--n_rep: The number of repetitions of predictions, each one with a unique cross-validation scheme<br>
--n_fold: The number of cross-validation folds<br>
--output: Which file to output the results to<br>
The following arguments are only important if the additional, reduced SNP-set prediction capability is used.<br>
--sigfile: The SNP-set file containing the SNP subsets.
--gwas_trait: The name of the SNP subset in the SNP-set file to use.
--n_random: The number of random SNP sets to use to predict for each n_rep.


## Example Usage

`Rscript Scripts/Predict.R --pheno 20260615_Lotus_ave_rep_family_iso_mds_rarefied_newformat.csv \<br>
                           --prediction_trait Seed_length                                       \<br>
                           --GD Lotus_GD.csv                                                    \<br>
                           --n_rep 100                                                          \<br>
                           --n_fold 5                                                           \<br>
                           --output gBLUP_Prediction_Seed_length.csv                            \<br>
                           --sigfile Empty_sigfile.csv                                          \<br>
                           --gwas_trait Seed_length                                             \<br>
                           --n_random 0`

## Outputs

A CSV file containing prediction accuracies, denoting the prediction type (full SNP set, reduced SNP set, random SNP set), the prediction accuracy, the repetition the accuracy belongs to, what the predicted trait is and what the name of SNP subset is.

## Runtime

Around an hour.
