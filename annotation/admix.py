###################
## set variables ##
###################

# specify file name prefix
in_name = 'admix'

# specify path to 1000g vcf files
ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'

# specify desired output directory
out_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/'

# set file paths to exectuables
plink_path = '/users/selasady/my_titan_itmi/tools/plink/plink --noweb'
tabix_path = '/users/selasady/my_titan_itmi/tools/tabix-0.2.6/tabix'
vcftools_path = '/titan/ITMI1/workspaces/users/selasady/tools/vcftools_0.1.13/bin/vcftools'
java_path = '/tools/java/jdk1.7/bin/java'
gatk_path = '/users/selasady/my_titan_itmi/tools/GenomeAnalysisTK.jar'

# import required modules
import os
import subprocess as sp
import pandas as pd

##############################
##  create known reference  ##
##############################
'''
downloaded 1000g filtered ped and map files from :
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/admixture_files/README.admixture_20141217
ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05.ped
ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05.map
filtered to retain biallelic, non-singleton SNV sites that are a minimum of 2KB apart from each other
merged the chromosomal files

convert ped to bed file and keep only maf between .05 and .5, with
LD pruning to inactivate first marker of pair within 50 marker window having CHM R^2 value > 0.5, steps of 5
see documentation http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
'''

# MAF and LD filtering
plink_file = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
plink_cmd = "{} --file {} --out all_chroms_100g_maf_ld --maf 0.05 --max-maf .5 --indep-pairwise 50 5 0.5".format(plink_path, plink_file)
# try:
#     ps = sp.Popen(plink_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#     print ps.communicate()[0]
# except sp.CalledProcessError as e:
#     print e.output

# retain only ld pruned variants and convert ped to bed/bim/fam file using plink
pruned_file = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_100g_maf_ld.prune.in'
ped2bed_cmd = '{} --file {} --extract {} --make-bed --out all_chroms_1000g_pruned'.format(plink_path, plink_file, pruned_file)
# try:
#     ps = sp.Popen(ped2bed_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#     print ps.communicate()[0]
# except sp.CalledProcessError as e:
#     print e.output

# convert the binary bed file to a readable map file to get genomic locations
#/users/selasady/my_titan_itmi/tools/plink/plink --noweb --bfile all_chroms_1000g_pruned --recode --out all_chroms_1000g_pruned_exported

## create .pop file

# read in .fam file and original 1000g ped file specifying ancestry
map_file = pd.read_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/integrated_call_samples_v2.20130502.ALL.ped', sep='\t')
fam_cols = ['family', 'subject', 'paternal', 'maternal', 'sex','affection']
fam_file = pd.read_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.fam', sep=' ', names = fam_cols)

# match .fam file with 1000g  integrated_call_samples_v2.20130502.ALL.ped by subject id
pop_df = fam_file.merge(map_file, left_on='subject', right_on='Individual ID', how = 'left')

print(pop_df.sample(20))


# # keep just the population column and save to file
# pop_out = pop_df['Population']
# pop_out.to_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.pop', header=False, index=False)

############################################
##  extract marker regions from vcf files ##
############################################

# for each subject vcf
    # extract regions that match marker vcf
    # append to marker file


####################
##  run admixture ##
####################