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

# # create function to run bash command with subprocess
def subprocess_cmd(command, cwd=os.getcwd()):
    '''
    Run programs in bash via subprocess
    :param command: command string as would be run on the command line
    :param input_dir: optional directory to run command in, default cwd
    :return: runs bash command
    '''
    print ("Running \n {}".format(command))
    ps = sp.Popen(command, shell=True,stdout=sp.PIPE,stderr=sp.PIPE, cwd=cwd)
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

def strip_vcf(input_dir):
    '''
    process all vcf.gz files to retain only chrom, pos, ref, alt and gt
    :param input_dir: path to directory where vcf files are located
    :return: _verified.vcf.gz files to merge with marker file before running ADMIXTURE
    '''
    for file in os.listdir(vcf_dir):
        if file.endswith('.vcf.gz'):
            print "Verifying VCF format for {}... \n".format(file)
            vcf_checked_out = str('.'.join(file.split('.')[:-2]) if '.' in file else file) + '_verified.vcf.gz'
            snp_verify_cmd = 'zcat {}/{} | {} | gzip > {}/{} '.format(vcf_dir, file, vcf_verify,vcf_dir, vcf_checked_out)
            ps = sp.Popen(snp_verify_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
            print ps.communicate()[0]

def vcf_to_bed(input_dir):
    '''
    convert stripped vcf to bed/bim/fam with plink 1.9
    :param input_dir: directory containing _verified.vcf.gz created from strip_vcf()
    :return:
    '''
    for file in os.listdir(input_dir):
        if file.endswith('_verified.vcf.gz'):
            vcf_out = str('.'.join(file.split('.')[:-2]) if '.' in file else file)
            vcf2plink_cmd = "{} --vcf {}/{} --double-id --biallelic-only strict --memory 300000 --geno 0.1 --allow-no-sex --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out {}/{}".format(plink_path, input_dir, file, vcf_dir, vcf_out)
            print ("Converting {} from vcf to bed/bim/fam").format(file)
            ps = sp.Popen(vcf2plink_cmd,shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
            print ps.communicate()[0]


def process_vcf():
    strip_vcf(vcf_dir)
    vcf_to_bed(vcf_dir)

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
def findBedStems(input_dir, file_ending):
    '''
    Locate bed files in a dir and ensure they have matching bim/fam
    :param input_dir: dir where bed files to merge are located
    :param file_ending: typicall '.bed'
    :return: list of bed files to merge
    '''
    bedList = []
    for file in os.listdir(input_dir):
        ref_bed = "{}.bed".format(ped_base)
        if (file.endswith(file_ending) or file.endswith(ref_bed)):
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
    out_file = "{}/{}".format(vcf_dir, timeStamped('_merge_files.txt'))
    with open(out_file,'wb') as outf:
        writer = csv.writer(outf, delimiter='\t')
        writer.writerows(file_list)

def plink2_merge(ref_dir):
    '''
    Generate command to merge bed files using plink then run with subprocess
    :param plink_path: path to plink exe
    :param marker_no_extension: location of marker file with no extension
    :param merge_file: path to merge text file create with create_merge_text()
    :param out_file: desired file name for output merged file
    :return: merged bed/bim/fam file of 1000g marker plus input bed files
    '''
    print ("Attempting to merge bed files with marker file.")
    merge_file = "{}/{}".format(vcf_dir, timeStamped('_merge_files.txt'))
    out_file = "{}/{}".format(vcf_dir, timeStamped('_merged'))
    merge_cmd = "{} --merge-list {} --make-bed --memory 40000 --out {}".format(plink_path, merge_file, out_file)
    subprocess_cmd(merge_cmd, vcf_dir)

def merge_vcf():
    '''
    run pipeline to merge vcf files with marker set
    :return: merged bed/bim/fam file to run through ADMIXTURE program
    '''
    #copy_marker(ref_dir, ped_base, vcf_dir)
    # merged_file = findBedStems(vcf_dir, '_verified.bed')
    # create_merge_text(merged_file)
    plink2_merge(ref_dir)


########################################################

if __name__ == '__main__':

    # path to input vcf files
    vcf_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2'
    # path to reference ped/map files
    ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf'
    # ped/map file basename
    ped_base = 'ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
    # panel file basename
    ref_panel = 'integrated_call_samples_v3.20130502.ALL.panel'

    # mark true if reference marker needs to be created, else False
    create_marker = 'True'

    # if create_marker == 'True':
    #     print ("Making marker... ")
    #     make_marker_function()

    #process_vcf()

    merge_vcf()

    #plink2_merge(ref_dir)




