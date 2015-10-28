# Impala Scripts
Scripts for accessing and querying impala. 

### [ACMG](https://github.com/summerela/impala_scripts/tree/master/ACMG) 
Script for locating variants in ACMG gene regions

###  [Annotation](https://github.com/summerela/impala_scripts/tree/master/annotation)
#### [admix.ipynb](https://github.com/summerela/impala_scripts/blob/master/annotation/admix.ipynb): script for running the ADMIXTURE program and uploading the results to impala
#### [global_distinct_variants.ipynb](https://github.com/summerela/impala_scripts/blob/master/annotation/global_distinct_variants.ipynb): 
Script to add the following annotations to variants in the global variants table:  
    - SnpEff coding consequence preditions  
    - DANN scores  
    - Ensembl gene annotations  
    - dbSNP rsid  
    - Clinvar clinical significance rating  
    - Kaviar allele frequency  
#### [snpeff.ipynb](https://github.com/summerela/impala_scripts/blob/master/annotation/snpeff.ipynb): run snpeff on each variant in the global variants table and add results as a table in impala
### [Clarity](https://github.com/summerela/impala_scripts/tree/master/clarity)
Scripts run on the clarity2 challenge data
### [ClinVar](https://github.com/summerela/impala_scripts/tree/master/ClinVar)
#### [clinvar_pathogenic.ipynb](https://github.com/summerela/impala_scripts/blob/master/ClinVar/clinvar_pathogenic.ipynb)
Script for identifying variants marked as pathogenic in ClinVar
### [NBS_genes](https://github.com/summerela/impala_scripts/tree/master/nbs_genes)
Contains versions of the script ran on data sets for locating variants in newborn screening gene regions
#### [nbs_pipeline](https://github.com/summerela/impala_scripts/blob/master/nbs_genes/nbs_pipeline..ipynb)
Script for locating variants in newborn screening gene regions
### [Parse_impala](https://github.com/summerela/impala_scripts/tree/master/parse_impala)
Helper functions for parsing impala and running queries
### [Scripts](https://github.com/summerela/impala_scripts/tree/master/Scripts)
Contains the following basic scripts
### [calc_maf.py](https://github.com/summerela/impala_scripts/blob/master/Scripts/calc_maf.py)
Calculate minor allele frequency in python
#### [find_illumina_vars.ipynb](https://github.com/summerela/impala_scripts/blob/master/Scripts/find_illumina_vars.ipynb)
A generic script for locating variants with an input list of samples and genes to look in
#### [noramlize.py](https://github.com/summerela/impala_scripts/blob/master/Scripts/normalize.py)
Perl script written by Terri Farrah to normalize variants that Joe Slagel will be converting to python
### [Tutorials](https://github.com/summerela/impala_scripts/tree/master/tutorials)
Tutorials for accessing and querying impala
