# Impala Scripts
Scripts for accessing and querying impala. 

## ACMG 
Script for locating variants in ACMG gene regions

##  Annotation
### admix.ipynb: script for running the ADMIXTURE program and uploading the results to impala
### global_distinct_variants.ipynb: script to add:  

    - SnpEff coding consequence preditions  
    - DANN scores  
    - Ensembl gene annotations  
    - dbSNP rsid  
    - Clinvar clinical significance rating  
    - Kaviar allele frequency  
### snpeff.ipynb: run snpeff on each variant in the global variants table and add results as a table in impala
## Clarity
Scripts run on the clarity2 challenge data
## ClinVar 
###clinvar_pathogenic.ipynb
Script for identifying variants marked as pathogenic in ClinVar
## NBS_genes
Contains versions of the script ran on data sets for locating variants in newborn screening gene regions
### nbs_pipeline
Script for locating variants in newborn screening gene regions
## Parse_impala
Helper functions for parsing impala and running queries
## Scripts
Contains the following basic scripts
### calc_maf.py 
Calculate minor allele frequency in python
### find_illumina_vars.ipynb
A generic script for locating variants with an input list of samples and genes to look in
### noramlize.py
Perl script written by Terri to normalize variants that Joe will be converting to python
## Tutorials
Tutorials for accessing and querying impala