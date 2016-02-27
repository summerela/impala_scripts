###################
## set variables ##
###################

# specify path to 1000g vcf files
#ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'
ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/testing'
kg_ref = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.bed'
# path to 1000g ped file
plink_file = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'

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
import csv
import datetime

##############################
##  create known reference  ##
##############################

### MAF and LD filtering ###



# command to create filter file 1000g ped file for MAF between .05 and .5, and LD 50,5,0.5
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

### create .pop file ###

# read in .fam file and original 1000g ped file specifying ancestry
# map_file = pd.read_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/integrated_call_samples_v2.20130502.ALL.ped', sep='\t')
# fam_cols = ['family', 'subject', 'paternal', 'maternal', 'sex','affection']
# fam_file = pd.read_csv('/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.fam', sep=' ', names = fam_cols)
#

# function to craete pop files from .fam and .ped
def make_pop(fam_file, map_file, out_pop):
    pop_df = fam_file.merge(map_file, left_on='subject', right_on='Individual ID', how = 'left')
    pop_out =  pop_df['Population']
    pop_out.to_csv(out_pop, header=False, index=False)

# make_pop(fam_file, map_file, "all_chroms_1000g_pruned.pop")

############################################
##  extract marker regions from vcf files ##
############################################

# for each subject vcf
    # extract regions that match marker vcf
    # TODO get code from Denise

# # convert filtered vcf to ped/map
# os.chdir(ref_dir)
# for file in os.listdir(ref_dir):
#     if file.endswith('.vcf.gz'):
#         vcf_out = str('.'.join(file.split('.')[:-2]) if '.' in file else file)
#         vcf2plink_cmd = "{} --gzvcf {} --plink --out {}".format(vcftools_path, file, vcf_out)
#         print ("Converting {} from vcf to ped/map").format(file)
#         try:
#             sp.call(vcf2plink_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.output
#
# # convert ped/map to bed/bim/fam
# for file in os.listdir(ref_dir):
#     if file.endswith('.ped'):
#         ped_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file)
#         bed_cmd = '{} --file {} --make-bed --allow-no-sex --out {}'.format(plink_path, ped_out, ped_out)
#         print ("Converting {} from ped/map to binary ped").format(ped_out)
#         try:
#             sp.call(bed_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.output

#########################
### merge bed files  ####
#########################
os.chdir(ref_dir)

# function to add time stamp to output file name
def timeStamped(fname, fmt='%Y-%m-%d{fname}'):
    return datetime.datetime.now().strftime(fmt).format(fname=fname)

# function locate bed files to merge
def findBedStems(input_dir, file_ending):
        """Find Plink binary stems (paths without .bed, .bim, .fam extension)"""
        bedList = []
        for file in os.listdir(input_dir):
            if (file.endswith(file_ending) and not file.startswith('all_chroms')):
                bedList.append(file)
        bedList.sort()
        plinkList = []
        for bed in bedList:
            (stem, ext) = os.path.splitext(bed)
            plinkList.append(stem)
            for suffix in (".bim", ".fam"):
                myFile = stem+suffix
                if not os.path.exists(input_dir + '/' + myFile):
                    print ("Missing Plink data file "+myFile)
        return plinkList

# function to create txt file of files to merge for plink input
def create_merge_text(merge_list):
    file_list = []
    for item in merge_list:
        bed_out = str(item + '.bed')
        bim_out = str(item + '.bim')
        fam_out = str(item + '.fam')
        t = [bed_out, bim_out, fam_out]
        file_list.append(t)
    with open(timeStamped('merge_files.txt'),'a') as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerows(file_list)

# define function to copy marker file to current wd
def copy_marker(marker_no_extension, current_path):
    copy_marker_cmd = "cp {}.* {}".format(marker_no_extension, current_path)
    try:
        sp.call(copy_marker_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    except sp.CalledProcessError as e:
        print e.output

kg_ref_base = str('.'.join(kg_ref.split('.')[:-1]) if '.' in kg_ref else kg_ref)
# copy_marker(kg_ref_base, ref_dir)

# function to merge all files listed in merge_list
def merge_bed(plink_path, marker_file, merge_list, out_file):
    merge_cmd = "{} --bfile {} --merge-list {} --make-bed --allow-no-sex --out {}".format(plink_path, marker_file, merge_list, out_file)
    try:
        sp.call(merge_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    except sp.CalledProcessError as e:
        print e.output

# function to remove mis-matching snps
def remove_snps(in_dir, plink_path, missnp_file):
     for file in os.listdir(ref_dir):
        base_name = str('.'.join(file.split('.')[:-1]) if '.' in file else file)
        try:
            remflip_cmd = "{} --bfile {} --exclude {} --make-bed --out {}_extracted".format(plink_path, base_name, missnp_file, base_name)
            print ("Removing mis-matched snps from {}").format(base_name)
            sp.call(remflip_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
        except sp.CalledProcessError as e:
            print e

# define function to remove file
def remove_file(file_to_remove):
    os.remove(file_to_remove)
    print("Deleted {} file.").format(file_to_remove)

# if there are erroneous snps, remove them and try the merge again
missnp_file = timeStamped('_merged') + '.missnp'





# else:
#     pass



## TODO test above, re-run merge on output, see if any additional errors need to be handled, add to process, rerun merge


# add samples to pop file






####################
##  run admixture ##
####################






###### TO DO #########
# TODO move reference and associcated index files to a common folder
# TODO move marker file and associated files to common folder
# TODO organize script




######## unused  ############
# dl ref.fa from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

## created input metadata file with following query
# select distinct family, genome, 'unknown' as fatherid, 'unknown' as motherid, gender, '-9' as phenotype
# from sgi_mapping
# where subject in ('F-101-371', 'F-101-346' ,'F-101-379', 'F-101-388', 'NB-101-527')

# os.chdir(ref_dir)

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
#
# # locate bed files to merge
# def findBedStems(input_dir):
#         """Find Plink binary stems (paths without .bed, .bim, .fam extension)"""
#         bedList = []
#         for file in os.listdir(input_dir):
#             if (file.endswith('.bed') and not file.startswith('all_chroms')):
#                 bedList.append(file)
#         bedList.sort()
#         plinkList = []
#         for bed in bedList:
#             (stem, ext) = os.path.splitext(bed)
#             plinkList.append(stem)
#             for suffix in (".bim", ".fam"):
#                 myFile = stem+suffix
#                 if not os.path.exists(myFile):
#                     print ("Missing Plink data file "+myFile)
#         return plinkList
#
# merge_list = findBedStems(ref_dir)
#
# # create txt file of files to merge for plink input
# test_list = []
#
# for item in merge_list:
#     bed_out = str(item + '.bed')
#     bim_out = str(item + '.bim')
#     fam_out = str(item + '.fam')
#     t = [bed_out, bim_out, fam_out]
#     test_list.append(t)
#
#
def timeStamped(fname, fmt='%Y-%m-%d{fname}'):
    return datetime.datetime.now().strftime(fmt).format(fname=fname)
# #
# import csv
# with open(timeStamped('merge_files.txt'),'a') as outf:
#     writer = csv.writer(outf, delimiter='\t')
#     writer.writerows(test_list)
#
# merge bed/bim/fam files
# bfile_basename = str('.'.join(kg_ref.split('.')[:-1]) if '.' in kg_ref else kg_ref)
# merge_cmd = "{} --bfile {} --merge-list {} --make-bed --out {}".format(plink_path, bfile_basename, timeStamped('merge_files.txt'), timeStamped('_merged'))
# try:
#     sp.call(merge_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
# except sp.CalledProcessError as e:
#     print e.output
