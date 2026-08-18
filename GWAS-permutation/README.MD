## GWAS permutation

This script is used for performing GWAS permutation in an adaptive way, that is, it stops prematurely if a sufficient, predetermined amount of permutations have a more extreme result than the unpermuted statistic. This is done to save computing power, especially when many permutations are done and some variants are not very significant. The script is made to run as part of a pipeline, as it only runs on a single SNP-trait association, and taking precomputed inputs such as the eigendecomposition of the genomic relationship matrix.

## Inputs

`20260818_reproducility_test_snpfile.csv` [File with SNP name, Chromosome, Position and genotypes]<br>
`20260818_reproducility_test_1.phenotypes.txt` [File with numeric phenotypic values]<br>
`20260818_reproducility_test_1.eigenU.txt` [File containing the eigenvectors of the genomic relationship matrix]<br>
`20260818_reproducility_test_1.eigenD.txt` [File containing the eigenvalues of the genomic relationship matrix]<br>
`20260818_reproducility_test_1.accessions.txt` [File containing the accession names of the phenotype file]<br>

## Arguments

Arguments are:<br>
    gemma_path (str): Path of the GEMMA executable.<br>
    phenotype (str): Path of the phenotype file.<br>
    snpfile (str): Path of the genotype file.<br>
    accessions (str): Path of the accessions file.<br>
    eigenvec (str): Path of the eigenvector file.<br>
    eigenval (str): Path of the eigenvalue file.<br>
    covariates (str): Path of the covariates file (Optional).<br>
    temp_dir (str): Temporary directory for temporary, intermediate files.<br>
    snp (str): The name of the SNP to be tested.<br>
    output (str): Path of the output file.<br>
    n_perm (list): Space separated list of integers specifying the number of permutations to perform.<br>
    perm_stop (int): Number of times the association statistic must exceed the original
        association statistic for the script to stop prematurely.

## Outputs

A CSV file containing the following columns:<br>
The original p-value for the association between the SNP and trait,<br>
the p-value for the association between the SNP and trait after permutation testing,<br>
the number of permutations performed,<br>
and the number of association statistics in the null distribution that exceed the original association statistic.

## Example use

```bash
python ./Scripts/Adaptive_GWAS_permutation.py                   \
    --accessions 20260818_reproducibility_test_1.accessions.txt \
    --phenotype 20260818_reproducibility_test_1.phenotypes.txt  \
    --eigenvec 20260818_reproducibility_test_1.eigenU.txt       \
    --eigenval 20260818_reproducibility_test_1.eigenD.txt       \
    --snpfile 20260818_reproducibility_test_snpfile.csv         \
    --temp_dir ./TMP_DIR                                        \
    --snp SNP_1_20972639                                        \
    --output 20260818_reproducibility_test_1.SNP_1_20972639.permutation.csv
```

## Runtime

Varies widely depending on architecture, sample size, and significance - varying from miliseconds to days.