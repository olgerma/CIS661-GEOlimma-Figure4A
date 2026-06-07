# CIS661-GEOlimma-Figure4A Replication
## Author
Maya Olger

## Project Description
This project replicates Figure 4A from the GEOlimma article by using the GSE8052 asthma dataset from the NCBI Gene Expression Omnibus (GEO). The purpose is to compare GEOlimma and Limma performance using ROC/AUC analysis on multiple random samples. 

## Files Included

### finalfigure4a.R
Main R script used for the replication.

### Figure4A_results.csv
Contains the AUC results and AUC improvements generated during the replication.

### Figure4A_full_replication.png
Final Figure 4A replication plot showing AUC improvement across subset trials.

## Required Packages

The following R packages are required:

- GEOquery
- limma
- hgu133plus2.db
- AnnotationDbi
- pROC
- ggplot2

## Additional Files Required

The following supplementary GEOlimma files must be downloaded from the original GEOlimma publication:

- GEOlimma_prob.rda
- 12859_2020_3932_MOESM8_ESM.r

These files should be placed in the project folder before running the script.

## Run Order

1. Install and load all required R packages.
2. Download the GEOlimma supplementary files.
3. Place the supplementary files in the project folder.
4. Open and run `finalfigure4a.R`.
5. Review the generated results and figure outputs.

## Reproducibility

A random seed (`set.seed(661)`) was used in order to make sure that the subset samples are reproducible. There are 40 trials for each subset size (n = 6, 9, 12, and 15) and calculates AUC improvement as:

AUC Improvement = GEOlimma AUC − Limma AUC

## Dataset

GSE8052: Asthma vs. Non-Asthma

Source:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8052

Platform:
GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)

## Repository Purpose

This repository is an implementation of the Figure 4A GEOlimma replication experiment completed for CIS 661. 
