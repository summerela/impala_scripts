###################
## set variables ##
###################
# TODO move under __main__
# test vcf files in /users/dmauldin/admixture

# specify path to 1000g vcf files
#vcf_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'
vcf_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2'
ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf'
kg_ref = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_1000g_pruned.bed'
# path to 1000g ped file
plink_file = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
ref_panel = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/integrated_call_samples_v3.20130502.ALL.panel'

# set file paths to exectuables
plink_path = '/users/selasady/my_titan_itmi/tools/plink1.9/plink --noweb'
vcftools_path = '/titan/ITMI1/workspaces/users/selasady/tools/vcftools_0.1.13/bin/vcftools'
java_path = '/tools/java/jdk1.7/bin/java'
gatk_path = '/users/selasady/my_titan_itmi/tools/GenomeAnalysisTK.jar'
admixture_path = '/users/selasady/my_titan_itmi/tools/admixture_linux-1.3.0/admixture'
vcf_verify = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/vcfBareBones.pl'

# import required modules
import os
import subprocess as sp
import pandas as pd
import csv
import datetime

##############################
##  create known reference  ##
##############################

# create function to run bash command with subprocess
# def run_bash_cmd(input_cmd):
#     ps = sp.Popen(input_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#     try:
#         ps.communicate()[0]
#     except sp.CalledProcessError as e:
#          print e.output



### create sort file to sort output by population ###

# change to dir where ref files are located
os.chdir(ref_dir)

# read in panel file and sort by population
pop_file = pd.read_csv(ref_panel, sep='\t')
pop_file.sort_values(by=['pop'], axis=0, inplace=True)
# pop_file.to_csv('sorted_ref.panel', sep='\t', header=False, index=False)

# subset panel file to contain only sample id's for sorting plink output
# pop_sort = pd.DataFrame()
# pop_sort['family_id'] = pop_file['sample']
# pop_sort['ind_id'] = pop_file['sample']
# pop_sort.to_csv('1000g_sort.txt', index=False, sep='\t', header=None)


### MAF and LD filtering ###

# change to directory with reference files

# command to create filter file 1000g ped file for MAF between .05 and .5, and LD 50,5,0.5
plink_cmd = "{} --file {} --out all_chroms_100g_maf_ld --maf 0.05 --max-maf .49 --indep-pairwise 50 5 0.5".format(plink_path, plink_file)
# retain only ld pruned variants and convert ped to bed/bim/fam file using plink
pruned_file = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/all_chroms_100g_maf_ld.prune.in'
ped2bed_cmd = '{} --file {} --extract {} --make-bed --indiv-sort file ./1000g_sort.txt --out all_chroms_1000g_pruned'.format(plink_path, plink_file, pruned_file)
# convert the binary bed file to a readable map file to get genomic locations for extracting vcf regions
map_pos_cmd = "/users/selasady/my_titan_itmi/tools/plink/plink --noweb --bfile all_chroms_1000g_pruned --recode --out all_chroms_1000g_pruned_exported"

# TODO don't exectue one subp until one passes
run_bash_cmd(plink_cmd)
run_bash_cmd(ped2bed_cmd)
run_bash_cmd(map_pos_cmd)

############################################
##  extract marker regions from vcf files ##
############################################

# for each subject vcf
    # extract regions that match marker vcf
    # TODO get code used from Denise


#############################
### Pre-process VCF Files ###
#############################
# process all vcf files to retain only chrom, pos, ref, alt and gt
# for file in os.listdir(vcf_dir):
#     if file.endswith('.vcf.gz'):
#         print "Verifying VCF format for {}... \n".format(file)
#         vcf_checked_out = str('.'.join(file.split('.')[:-2]) if '.' in file else file) + '_verified.vcf'
#         snp_verify_cmd = 'zcat {} | {} | gzip -c > {}.gz '.format(file, vcf_verify,vcf_checked_out)
#         # create the file and run snpeff
#         with open(vcf_checked_out, "w") as out_file:
#             try:
#                 ps = sp.Popen(snp_verify_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#                 print ps.communicate()[0]
#             except sp.CalledProcessError as e:
#                  print e.output

# TODO above pipe is creating empty temp file, fix

# convert vcf to bed/bim/fam with plink 1.9
# for file in os.listdir(vcf_dir):
#     if file.endswith('_verified.vcf.gz'):
#         vcf_out = str('.'.join(file.split('.')[:-2]) if '.' in file else file)
#         vcf2plink_cmd = "{} --vcf {} --double-id --biallelic-only strict --geno 0.1 --allow-no-sex --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out {}".format(plink_path, file, vcf_out)
#         print ("Converting {} from vcf to bed/bim/fam").format(file)
#         try:
#             sp.call(vcf2plink_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.message


######################################
### functions to merge bed files  ####
######################################

# define function to copy marker file to current wd
def copy_marker(marker_no_extension, current_path):
    '''
    copy 1000g bed/bim/fam marker files to same dir as vcf files to merge
    :param marker_no_extension: name of marker file without file extension
    :param current_path: vcf_dir to copy marker file to
    :return: copies marker bed/bim/fam file to directory where vcf files are located
    '''
    copy_marker_cmd = "cp {}.* {}".format(marker_no_extension, current_path)
    try:
        sp.call(copy_marker_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    except sp.CalledProcessError as e:
        print e.output

# add time stamp to output file name
def timeStamped(fname, fmt='%Y-%m-%d{fname}'):
    '''
    add a time stamp to output files
    :param fname: file name to stamp
    :param fmt: time stamp format
    :return: file named with todays date
    '''
    return datetime.datetime.now().strftime(fmt).format(fname=fname)

# function locate bed files to merge
def findBedStems(input_dir, file_ending):
    '''
    Locate bed files in a dir and ensure they have matching bim/fam
    :param input_dir: dir where bed files to merge are located
    :param file_ending: typicall '.bed'
    :return: list of bed files to merge
    '''
    bedList = []
    for file in os.listdir(input_dir):
        ref_name = os.path.basename(kg_ref)
        if (file.endswith(file_ending) or file == ref_name):
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
    '''
    create a three column text file of bed/bim/fam files to merge based
    on list from findBedStems
    :param merge_list: list generated with findBedStems()
    :return: tsv file of each bed bim fam file to  merge for input to plink merge
    '''
    file_list = []
    for item in merge_list:
        bed_out = str(item + '.bed')
        bim_out = str(item + '.bim')
        fam_out = str(item + '.fam')
        t = [bed_out, bim_out, fam_out]
        file_list.append(t)
    with open(timeStamped('_merge_files.txt'),'wb') as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerows(file_list)

def plink2_merge(plink_path, merge_file, out_file):
    '''
    Generate command to merge bed files using plink then run with subprocess
    :param plink_path: path to plink exe
    :param marker_no_extension: location of marker file with no extension
    :param merge_file: path to merge text file create with create_merge_text()
    :param out_file: desired file name for output merged file
    :return: merged bed/bim/fam file of 1000g marker plus input bed files
    '''
    print ("Attempting to merge bed files with marker file.")
    merge_cmd = "{} --merge-list {} --make-bed --out {}".format(plink_path, merge_file, out_file)
    try:
        sp.call(merge_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    except sp.CalledProcessError as e:
        print e.output



# define function to remove file
def remove_file(file_to_remove):
    '''
    delete missnp file so that snp removal function stops running once
    all missing snips have been removed
    :param file_to_remove: path to missnp file
    :return: directory with no missnp file
    '''
    os.remove(file_to_remove)
    print("Deleted {} file.").format(file_to_remove)

# function to remove mis-matching snps
def remove_snps(in_dir, plink_path, missnp_file):
    '''
    remove snps listed in mssnp file from input bed files before
    trying to remerge
    :param in_dir: directory where input bed files are located
    :param plink_path: path to plink executable
    :param file_no_extenstion: path to marker file with no file extension
    :param missnp_file: path to missnp file
    :return: bed files with missing snps removed
    '''
    for file in os.listdir(in_dir):
        if file.endswith('.bed'):
            base_name = str('.'.join(file.split('.')[:-1]) if '.' in file else file)
            try:
                remflip_cmd = "{} --bfile {} --exclude {} --make-bed --out {}".format(plink_path, base_name, missnp_file, base_name)
                print ("Removing mis-matched snps from {}").format(base_name)
                sp.call(remflip_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
            except sp.CalledProcessError as e:
                print e

#####################
## merge bed files ##
#####################

# create base file name
kg_ref_base = str('.'.join(kg_ref.split('.')[:-1]) if '.' in kg_ref else kg_ref)

# # copy 1000g marker bed/bim/fam to vcf_dir
# copy_marker(kg_ref_base, vcf_dir)
#
#
# # find bed files to merge
# first_merge = findBedStems(vcf_dir, '_verified.bed')
# # create time-stamped merge text
# create_merge_text(first_merge)
#
# # attempt to merge files
# plink2_merge(plink_path, timeStamped('_merge_files.txt'), timeStamped('_merged'))
#
#
# # if there are erroneous snps, remove them and try the merge again
# missnp_file = timeStamped('_merged-merge') + '.missnp'
#
# # # while there is an error file
# while os.path.isfile(missnp_file):
#     # remove mis-matched snps from files
#     print ("Removing mis-matching SNP's... \n")
#     remove_snps(vcf_dir, plink_path, missnp_file)
#     # delete the file
#     remove_file(missnp_file)
#     # try merging again
#     plink2_merge(plink_path, timeStamped('_merge_files.txt'),timeStamped('_merged'))

def make_pop(ref_file, out_pop, kval):
    '''
    Create .pop file with known populations from 1000g and '-' for unknown from vcf files
    :param ref_file: .ped file from marker reference
    :param out_pop: base name for output .pop file
    :param kval: int of 4 or 26 to set kvalue for running admixture and creating pop file
    :return .pop file for merged bed/bim/fam set needed for running supervised admixture
    '''
    # pull population from the 1000g map file
    ref_map = pd.read_csv(ref_file, sep='\t')
    # subset for needed columns
    map_match =  ref_map[['Individual ID', 'Population']]
    # read in -merge.fam file
    for file in os.listdir(vcf_dir):
        if file.endswith('_merged.fam'):
            merged_fam = pd.read_csv(file, sep=' ', names = ['fam_id', 'Individual ID', 'Father', 'Mother', 'Gender', 'Phenotype'])
            merged_df = merged_fam[['Individual ID']]
            out_df = pd.merge(merged_df, map_match, how='left', on= 'Individual ID')
            out_df = out_df['Population']
            out_df.to_csv(out_pop, sep='\t', na_rep='-', header=False, index=False)

#make_pop(ref_panel, timeStamped('_merged.pop'))

####################
##  run admixture ##
####################

# create admix command with -j20 = 20 cores and --supervised 26 = 26 populations
admix26_cmd = "nohup {} -j20 {} --supervised 26".format(admixture_path, timeStamped("_merged.bed"))


# TODO subp is not calling the process, ran manually
# print admix26_cmd
# run_bash_cmd(admix26_cmd)

###############################
### upload results to impala ##
###############################



#############################
######## unused  ############
#############################
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


# convert vcf to ped/map
# for file in os.listdir(vcf_dir):
#     if file.endswith('.vcf.gz'):
#         vcf_out = str('.'.join(file.split('.')[:-2]) if '.' in file else file)
#         vcf2plink_cmd = "{} --gzvcf {} --plink --biallelic-only strict --out {}".format(vcftools_path, file, vcf_out)
#         print ("Converting {} from vcf to ped/map").format(file)
#         try:
#             sp.call(vcf2plink_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.message

# convert ped/map to bed/bim/fam
# for file in os.listdir(vcf_dir):
#     if file.endswith('.ped'):
#         ped_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file)
#         bed_cmd = '{} --file {} --make-bed --allow-no-sex --out {}'.format(plink_path, ped_out, ped_out)
#         print ("Converting {} from ped/map to binary ped").format(ped_out)
#         try:
#             sp.call(bed_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#         except sp.CalledProcessError as e:
#             print e.message

# # function to merge all files listed in merge_list
# def merge_bed(plink_path, marker_no_extension, merge_file, out_file):
#     '''
#     Generate command to merge bed files using plink then run with subprocess
#     :param plink_path: path to plink exe
#     :param marker_no_extension: location of marker file with no extension
#     :param merge_file: path to merge text file create with create_merge_text()
#     :param out_file: desired file name for output merged file
#     :return: merged bed/bim/fam file of 1000g marker plus input bed files
#     '''
#     print ("Attempting to merge bed files with marker file.")
#     merge_cmd = "{} --bfile {} --merge-list {} --make-bed --allow-no-sex --out {}".format(plink_path, marker_no_extension, merge_file, out_file)
#     try:
#         sp.call(merge_cmd, shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
#     except sp.CalledProcessError as e:
#         print e.output