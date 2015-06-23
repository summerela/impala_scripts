## clinvary_query.py ##

Description
************
This script queries the impala system to find either cgi or illumina variants that fall within regions marked with
4 or 5 clinical significance in clinvar, but not those that are ever marked 2 or 3.

Input
*****
The following command line options:

    --chr       =  Enter chromosomes of interest as comma-separated list, no spaces, leave blank for all
    --zygosity  =  Enter zygosity of interest as comma-separated list, no spaces, leave blank for all
                   Choices= hom-alt,hom-ref,het,hom
    --genes     =  Enter genes of interest as comma-separated list, no spaces, leave blank for all
    --member    =  Enter member types of interest as comma-separated list, no spaces, leave blank for all
                   Choices = M,F,NB
    --sample_id =  Enter sample id's of interest as comma-separated list, no spaces, leave blank for all
                   NOTE: enter in format compatible with desired platform
    --platform  =  Enter either cgi or illumina

    Example:
    clinvary_query.py --chr=7,7 --zygosity=hom --genes=HOXB1,BRCA1,BRCA2 --member=M,F --sample_id=101-5456-NB \
    --platform=illumina

 Output
 *******
 A csv file will be saved to the current working directory titled "clinvar_results.csv".

 Help
 ****
 For any questions about running this script, please contact Summer at selasady@systemsbiology.org.