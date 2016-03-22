'''
Run ADMIXTURE
- create reference marker file
- convert VCF of interest to plink bed format
- run admixture
- output results to impala table

input =
    ref ped/map files
    ref panel: file (three columns, tsv, (subject | population | super-pop ))
    vcf_dir: all vcf files in this directory will be merged with marker file
output = admixture results... TODO add more info here
'''

# import modules
import os
import subprocess as sp
import shlex
import pandas as pd
import csv
import datetime

# program paths
tool_path = '/users/selasady/my_titan_itmi/tools/'
plink_path = "{}plink1.9/plink".format(tool_path)
java_path = '/tools/java/jdk1.7/bin/java -jar -Xmx16g'
gatk_path = "{}GenomeAnalysisTK.jar".format(tool_path)
admixture_path = "{}admixture_linux-1.3.0/admixture".format(tool_path)
vcf_verify = "{}snpEff/scripts/vcfBareBones.pl".format(tool_path)


####################################
### Create reference marker file ###
####################################

# create text file with subject id's grouped by population for sorting
def create_sorter(ref_panel):
    '''
    :param ref_panel: panel file (three columns, tsv, col headers = sample | pop | super_pop )
    :return sorted_ref.txt a two column tsv ( sample | pop ) to be fed to plink for sorting by pop
    '''
    panel_file = "{}/{}".format(ref_dir, ref_panel)
    out_panel = "{}/sorted_panel.txt".format(ref_dir)
    out_sort = "{}/sorted_ref.txt".format(ref_dir)
    print("Saving sorted panel file as {}".format(out_panel))
    panel = pd.read_csv(panel_file, sep='\t')
    panel_out = panel[['sample', 'pop', 'super_pop', 'gender']]
    panel_out.sort_values(by=['pop'], axis = 0, inplace=True)
    panel_out.to_csv(out_panel, index=False, sep='\t')
    sort_out = panel_out[['sample', 'sample']]
    print ("Saving sorted reference text as {}".format(out_sort))
    sort_out.to_csv(out_sort, index=False, sep='\t', header=False)

# create function to run bash command with subprocess
def subprocess_cmd(command, input_dir):
    '''
    Run programs in bash via subprocess
    :param command: command string as would be run on the command line
    :param input_dir: optional directory to run command in, default cwd
    :return: runs bash command
    '''
    print ("Running \n {}".format(command))
    ps = sp.Popen(command, shell=True,stdout=sp.PIPE,stderr=sp.PIPE, cwd=input_dir)
    try:
       print ps.communicate()
    except sp.CalledProcessError as e:
         print e

def make_marker_function():
    '''
    Runs commands to create marker file if make_marker == True
    :return: marker file to merge with input VCF files for running ADMIXTURE
    '''
     # create sorted panel file
    create_sorter(ref_panel)
    # filter reference ped file to retain MAF between .05 and 0.5 and LD 50 5 0.5
    plink_cmd = "{} --file {}/{} --out {}/{}_maf_ld --maf 0.05 --max-maf .49 --indep-pairwise 50 5 0.5".format(plink_path, ref_dir, ped_base, ref_dir, ped_base)
    # retain only maf/ld pruned variants and convert ped to bed/bim/fam file using plink
    pruned_file = '{}/{}_maf_ld.prune.in'.format(ref_dir, ped_base)
    ped2bed_cmd = '{} --file {}/{} --extract {} --make-bed --indiv-sort file {}/sorted_ref.txt --out {}/{}'.format(plink_path, ref_dir, ped_base, pruned_file, ref_dir, ref_dir, ped_base)
    # convert the binary bed file to a readable map file to get genomic locations for extracting vcf regions
    map_pos_cmd = "{} --noweb --bfile {}/{} --recode --out {}/{}_map".format(plink_path, ref_dir, ped_base, ref_dir, ped_base)
    marker_cmds = [plink_cmd,ped2bed_cmd, map_pos_cmd]
    for cmd in marker_cmds:
        subprocess_cmd(cmd, ref_dir)

###############################
### Pre-process VCF file(s) ###
###############################

def strip_vcf(input_dir, input_vcf):
    '''
    process all vcf.gz files to retain only chrom, pos, ref, alt and gt
    :param input_dir: path to directory where vcf files are located
    :param input_vcf: vcf file to process
    :return: _verified.vcf.gz files to merge with marker file before running ADMIXTURE
    '''
    print "Verifying VCF format for {}... \n".format(input_vcf)
    vcf_checked_out = str('.'.join(input_vcf.split('.')[:-2]) if '.' in input_vcf else input_vcf) + '_verified.vcf.gz'
    snp_verify_cmd = 'zcat {}/{} | {} | gzip > {}/{} '.format(input_dir, input_vcf, vcf_verify,input_dir, vcf_checked_out)
    ps = sp.Popen(snp_verify_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    print ps.communicate()[0]

def vcf_to_bed(input_dir, input_vcf):
    '''
    convert stripped vcf to bed/bim/fam with plink 1.9
    :param input_dir: directory containing _verified.vcf.gz created from strip_vcf()
    :param input_vcf: vcf file to process
    :return: stripped, verified vcf files converted to plink binary (bed) format
    '''
    vcf_out = str('.'.join(input_vcf.split('.')[:-2]) if '.' in input_vcf else input_vcf)
    vcf2plink_cmd = "{} --vcf {}/{} --double-id --biallelic-only strict --memory 300000 --geno 0.1 --allow-no-sex --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out {}/{}".format(plink_path, input_dir, input_vcf, input_dir, vcf_out)
    print ("Converting {} from vcf to bed/bim/fam").format(input_vcf)
    ps = sp.Popen(vcf2plink_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    print ps.communicate()[0]


def process_vcf(input_dir):
    for file in os.listdir(input_dir):
                if file.endswith('vcf.gz'):
                    strip_vcf(input_dir, file)
                    if file.endswith('_verified.vcf.gz'):
                        vcf_to_bed(input_dir, file)

#########################
### merge bed files  ####
#########################

# define function to copy marker file to current wd
def copy_marker(ref_dir, marker_no_extension, vcf_dir):
    '''
    copy 1000g bed/bim/fam marker files to same dir as vcf files to merge
    :param marker_no_extension: name of marker file without file extension
    :param current_path: vcf_dir to copy marker file to
    :return: copies marker bed/bim/fam file to directory where vcf files are located
    '''
    print "Copying markfer file to {}".format(vcf_dir)
    copy_marker_cmd = "cp {}/{}.* {}".format(ref_dir, marker_no_extension, vcf_dir)
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
def findBedStems(input_dir, input_file, ped_base):
    '''
    Locate bed files in a dir and ensure they have matching bim/fam
    :param input_dir: dir where bed files to merge are located
    :param input_file: file to merge with marker file
    :param ped_base: basename of reference marker file to merge with vcf
    :return: list of bed files to merge
    '''
    bedList = []
    ref_bed = "{}.bed".format(ped_base)
    for file in os.listdir(input_dir):
        if file == ref_bed:
             bedList.append(file)
    bedList.append(input_file)
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

def create_merge_text(input_vcf, merge_list):
    '''
    create a three column text file of bed/bim/fam files to merge based
    on list from findBedStems
    :param input_vcf: vcf file to merge with marker file
    :param merge_list: list generated with findBedStems()
    :return: tsv file of each bed bim fam file to  merge for input to plink merge
    '''
    print ("Creating merge file for {}".format(input_vcf))
    file_list = []
    for item in merge_list:
        bed_out = str(item + '.bed')
        bim_out = str(item + '.bim')
        fam_out = str(item + '.fam')
        t = [bed_out, bim_out, fam_out]
        file_list.append(t)
    vcf_base = str('.'.join(input_vcf.split('.')[:-2]) if '.' in input_vcf else input_vcf)
    out_file = "{}/{}".format(vcf_dir, timeStamped('_{}_merge.txt'.format(vcf_base)))
    with open(out_file,'wb') as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerows(file_list)

def plink2_merge(input_vcf, vcf_dir):
    '''
    Generate command to merge bed files using plink then run with subprocess
    :param plink_path: path to plink exe
    :param marker_no_extension: location of marker file with no extension
    :param merge_file: path to merge text file create with create_merge_text()
    :param out_file: desired file name for output merged file
    :return: merged bed/bim/fam file of 1000g marker plus input bed files
    '''
    print ("Attempting to merge bed file with marker file.")
    vcf_base = str('.'.join(input_vcf.split('.')[:-2]) if '.' in input_vcf else input_vcf)
    merge_file = "{}/{}".format(vcf_dir, timeStamped('_{}_merge.txt'.format(vcf_base)))
    out_file = "{}/{}".format(vcf_dir, timeStamped('_{}_merged'.format(vcf_base)))
    merge_cmd = "{} --merge-list {} --make-bed --memory 40000 --out {}".format(plink_path, merge_file, out_file)
    subprocess_cmd(merge_cmd, vcf_dir)

def merge_vcf(vcf_dir):
    '''
    run pipeline to merge bed files with marker set
    :param vcf_dir: directory path to vcf files to merge
    :return: merged bed/bim/fam file to run through ADMIXTURE program
    '''
    # copy marker file to vcf_dir
    copy_marker(ref_dir, ped_base, vcf_dir)
    for file in os.listdir(vcf_dir):
        if file.endswith('_verified.bed'):
            merged_file = findBedStems(vcf_dir, file, ped_base)
            create_merge_text(file, merged_file)
            plink2_merge(file, vcf_dir)

##############################
### Check for merge errors ###
##############################

# function to remove mis-matching snps
def remove_snps(vcf_dir, input_vcf, plink_path, missnp_file):
    '''
    remove snps listed in mssnp file from input bed files before
    trying to remerge
    :param vcf_dir: directory where input bed files are located
    :param plink_path: path to plink executable
    :param missnp_file: path to missnp file
    :return: bed files with missing snps removed
    '''
    base_name = str('.'.join(input_vcf.split('.')[:-1]) if '.' in input_vcf else input_vcf)
    remflip_cmd = "{} --bfile {} --exclude {} --make-bed --out {}".format(plink_path, base_name, missnp_file, base_name)
    subprocess_cmd(remflip_cmd, vcf_dir)

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

def check_merge_errors(vcf_dir, input_vcf):
    # if there are erroneous snps, remove them and try the merge again
    vcf_base = str('.'.join(input_vcf.split('.')[:-2]) if '.' in input_vcf else input_vcf)
    missnp_file = "{}/{}.missnp".format(vcf_dir, timeStamped('_{}_merged'.format(vcf_base)))
    # while there is an error file
    while os.path.isfile(missnp_file):
        # remove mis-matched snps from files
        print ("Removing mis-matching SNP's... \n")
        remove_snps(vcf_dir, plink_path, missnp_file)
        # delete the file
        remove_file(missnp_file)
        # try merging again
        plink2_merge(plink_path, timeStamped('_merge_files.txt'),timeStamped('_merged'))

#######################
###  run admixture  ###
#######################

# make pop file with 5 superpops
def make_pop5(ref_file, out_pop, kval):
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

# make pop file with 26 populations



########################################################

if __name__ == '__main__':

    # path to input vcf files
    vcf_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test3'
    # path to reference ped/map files
    ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf'
    # ped/map file basename
    ped_base = 'ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
    # panel file basename
    ref_panel = 'integrated_call_samples_v3.20130502.ALL.panel'

    # mark true if reference marker needs to be created, else False
    create_marker = 'False'

    # if create_marker == 'True':
    #     print ("Making marker... ")
    #     make_marker_function()

    #process_vcf(vcf_dir)

    merge_vcf(vcf_dir)



# TODO add program path as arg to each function
# TODO consistent variable names as input args to functions
# TODO global var with vcf base name and pass to functions

