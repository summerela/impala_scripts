###################
## set variables ##
###################

# specify file name prefix
in_name = 'admix'

# specify path to 1000g vcf files
#ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'
ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/testing'

#GATK variables
# path to reference fasta file for gatk
ref_fasta = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/testing/hs37d5.fa'
metadata_file = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/testing/test_metadata.fam'

# specify desired output directory
out_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/'

# set file paths to exectuables
plink_path = '/users/selasady/my_titan_itmi/tools/plink/plink --noweb'
tabix_path = '/tools/bin/tabix'
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
# map_file = pd.read_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/integrated_call_samples_v2.20130502.ALL.ped', sep='\t')
# fam_cols = ['family', 'subject', 'paternal', 'maternal', 'sex','affection']
# fam_file = pd.read_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.fam', sep=' ', names = fam_cols)
#
# # match .fam file with 1000g  integrated_call_samples_v2.20130502.ALL.ped by subject id
# pop_df = fam_file.merge(map_file, left_on='subject', right_on='Individual ID', how = 'left')
#
# # # keep just the population column and save to file
# pop_out = pop_df['Population']
# pop_out.to_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.pop', header=False, index=False)

############################################
##  extract marker regions from vcf files ##
############################################

# for each subject vcf
    # extract regions that match marker vcf

# remove 'chr' prefix from chromosome column
# os.chdir(ref_dir)
# for i in *.vcf.gz;
#     do gunzip -c $i | perl -p -e 's/^chr//' | gzip -c > ${i%.vcf.gz}_stripped.vcf.gz;
# done

# dl ref.fa from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

## created input metadata file with following query
# select distinct family, genome, 'unknown' as fatherid, 'unknown' as motherid, gender, '-9' as phenotype
# from sgi_mapping
# where subject in ('F-101-371', 'F-101-346' ,'F-101-379', 'F-101-388', 'NB-101-527')

os.chdir(ref_dir)

# index the vcf files
# for file in os.listdir(os.getcwd()):
#     if file.endswith('.vcf.gz'):
#         tabix_cmd = "{} -p vcf {}".format(tabix_path, file)
#         print ("Running {}".format(tabix_cmd)
#         try:
#            sp.call(tabix_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.output

# convert stripped.vcf.gz to BED/BIM/FAM files

bedlist = []
bimlist= []
famlist = []

# for file in os.listdir(os.getcwd()):
#     if file.endswith('.vcf.gz'):
#         base_out = str('.'.join(file.split('.')[:-3]) if '.' in file else file)
#         bed_out = base_out + '.bed'
#         bim_out = base_out + '.bim'
#         fam_out = base_out + '.fam'
#         vcf2bed_cmd = str("{} -jar -Xmx16g {} -T VariantsToBinaryPed -R {} -V {} -m {} -bed {} -bim {} -fam {} -mgq 0"\
#                           .format(java_path, gatk_path, ref_fasta, file, metadata_file, bed_out, bim_out, fam_out))
#         print ("Running {}").format(vcf2bed_cmd)
#         try:
#             sp.call(vcf2bed_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.output

# copy marker file to current wd
copy_marker = "cp /users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.* {}".format(ref_dir)
# try:
#     sp.call(copy_marker, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
# except sp.CalledProcessError as e:
#     print e.output

# merge bed files
import datetime

# locate bed files to merge
def findBedStems(input_dir):
        """Find Plink binary stems (paths without .bed, .bim, .fam extension)"""
        bedList = []
        for file in os.listdir(input_dir):
            if (file.endswith('.bed') and not file.startswith('all_chroms')):
                bedList.append(file)
        bedList.sort()
        plinkList = []
        for bed in bedList:
            (stem, ext) = os.path.splitext(bed)
            plinkList.append(stem)
            for suffix in (".bim", ".fam"):
                myFile = stem+suffix
                if not os.path.exists(myFile):
                    print ("Missing Plink data file "+myFile)
        return plinkList

merge_list = findBedStems(ref_dir)

# create txt file of files to merge for plink input
test_list = []

for item in merge_list:
    bed_out = str(item + '.bed')
    bim_out = str(item + '.bim')
    fam_out = str(item + '.fam')
    t = [bed_out, bim_out, fam_out]
    test_list.append(t)


def timeStamped(fname, fmt='%Y-%m-%d{fname}'):
    return datetime.datetime.now().strftime(fmt).format(fname=fname)
#
import csv
with open(timeStamped('merge_files.txt'),'a') as outf:
    writer = csv.writer(outf, delimiter='\t')
    writer.writerows(test_list)

# merge bed/bim/fam files
merge_cmd = "{} --bfile all_chroms_1000g_pruned --merge-list {} --make-bed --out {}".format(plink_path, timeStamped('merge_files.txt'), timeStamped('_merged'))
try:
    sp.call(merge_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
except sp.CalledProcessError as e:
    print e.output


# add samples to pop file






# ./plink --file ALL.wgs.phase3_shapeit2_filtered.20141217.ped --out ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05 --make-bed

####################
##  run admixture ##
####################






###### TO DO #########
# TODO move reference and associcated index files to a common folder
# TODO move marker file and associated files to common folder
# TODO organize script