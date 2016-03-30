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
bcftools = "{}bcftools.bcftools".format(tool_path)
extract_script = "/titan/ITMI1/workspaces/users/dmauldin/extract_regions_stream_job.pl"


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

def make_marker(ref_panel, plink_path, ref_dir, ped_base):
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
    # # convert the binary bed file to a readable map file to get genomic locations for extracting vcf regions
    # map_pos_cmd = "{} --noweb --bfile {}/{} --recode 12 --out {}/{}_map".format(plink_path, ref_dir, ped_base, ref_dir, ped_base)
    # convert bed file to vcf
    bed_to_vcf_cmd = "{} --file {}  --recode vcf --out {}".format(plink_path, ped_base, ped_base)
    # extract chrom pos ref from marker vcf
    make_region = r'''cat {}.vcf  | grep -v '^#' | awk '{{OFS="\t";print "chr"$1, $2, $4}}' | gzip > {}_markers.gz'''.format(ped_base, ped_base)
    marker_cmds = [plink_cmd,ped2bed_cmd, bed_to_vcf_cmd, make_region]
    for cmd in marker_cmds:
        subprocess_cmd(cmd, ref_dir)

###############################
### Pre-process VCF file(s) ###
###############################

# subset metadata for input vcf file
def subset_metadata(metadata_file, input_vcf):
    '''
    subset metadata file specifying path and subject info for feeding to perl extract
    script; columns to include are the same as meta_names variable
    :param metadata_file: path to complete metadata file for each platform
    :param input_vcf: vcf file to subset. assumes vcf file begins 102-00029-03.LP6005637-DNA_H05.etc...
    as output by perl extract script
    :return: metadata subset dataframe containing vcf file info to feed to update_meta_path()
    '''
    vcf_basename = str('.'.join(input_vcf.split('.')[:1]) if '.' in input_vcf else input_vcf)
    #meta_names = ['Vendor',  'Bucket',  'Study',   'Family',  'Genome',  'Subject', 'Sample',  'Assembly', 'Gestalt ID', 'Gender',  'Term'
                  # 'Country of Birth', 'Member',  'VCF', 'PATH', 'S3','BIGDATA', 'GESTALT',  'ITMI_MAPPINGS', 'COUNT', 'S3URL', 'BigdataURL',      'GestaltURL', 'MappingsURL']
    # read in metadata file
    meta_file = pd.read_csv(metadata_file, sep='\t', header=0)
    # subset metadata file and update vcf file paths
    vcf_subset = meta_file[(meta_file['Genome'].str.contains(vcf_basename))]
    return vcf_subset

def update_meta_path(vcf_subset, input_dir, input_vcf, out_dir):
    '''
    Update vcf file path for subset metadata file for current vcf location
    :param vcf_subset: dataframe subset of platform metadata file created from subset_metadata()
    :param input_dir: directory containing vcf file
    :param input_vcf: vcf file to subset metadata for and update the file path to location of this vcf
    :param out_dir: directory to store output metadata text
    :return: subset metadata file containing information for input_vcf as tsv
    '''
    vcf_path = os.path.abspath(os.path.join(input_dir, input_vcf))
    print ("\n Creating metadata file for {}".format(vcf_path))
    vcf_subset['VCF'] = vcf_path
    vcf_basename = str('.'.join(input_vcf.split('.')[:1]) if '.' in input_vcf else input_vcf)
    out_file = "{}/{}_metadata.txt".format(out_dir, vcf_basename)
    vcf_subset.to_csv(out_file, sep='\t', index=False, na_rep="NA")

# extract regions from vcf file that match marker regions
def extract_variants(input_vcf, ped_base, ref_dir, out_dir):
    '''
    Runs extract_regions_stream_job.pl written by Denise Maulden to extract regions from input vcf files
     that match regions in the marker file to save computational time
    :param input_vcf: path to vcf file to extract regions from
    :param ped_base: base name of input ped file used to create marker file
    :param ref_dir: dictory path to reference marker file
    :param metadata_file: file with metadata and location of input vcf files made from create_metadata_file()
    :param out_dir: directory to write output vcf's
    :return: *_filtered.vcf.gz files retaining only regions that match the maker file regions
    '''
    print ("Extracting marker regions from vcf files... \n")
    vcf_basename = str('.'.join(input_vcf.split('.')[:1]) if '.' in input_vcf else input_vcf)
    region_file = "{}/{}_markers.gz".format(ref_dir, ped_base)
    metadata_file = "{}{}_metadata.txt".format(out_dir, vcf_basename)
    extract_cmd = "/tools/bin/perl {} --regionFile {} --metadata {} --outDir {}  --compressionType bzip2".format(extract_script, region_file, metadata_file, out_dir)
    subprocess_cmd(extract_cmd, vcf_dir)

def strip_vcf(input_dir, input_vcf):
    '''
    process all vcf.gz files to retain only chrom, pos, ref, alt and gt
    :param input_dir: path to directory where vcf files are located
    :param input_vcf: vcf file to process
    :return: _verified.vcf.gz files to merge with marker file before running ADMIXTURE
    '''
    print ("Verifying VCF format for {}... \n".format(input_vcf))
    vcf_checked_out = str('.'.join(input_vcf.split('.')[:-2]) if '.' in input_vcf else input_vcf) + '_verified.vcf.gz'
    snp_verify_cmd = 'bzcat {}/{} | {} | gzip > {}/{} '.format(input_dir, input_vcf, vcf_verify,input_dir, vcf_checked_out)
    subprocess_cmd(snp_verify_cmd, input_dir)

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
    subprocess_cmd(vcf2plink_cmd, input_dir)

# run process vcf functions
def extract_vars(metadata_file, input_dir, ref_dir, out_dir):
    for file in os.listdir(input_dir):
        if file.endswith('vcf.gz'):
            meta_subset = subset_metadata(metadata_file, file)
            update_meta_path(meta_subset, input_dir, file, out_dir)
            extract_variants(file, ped_base, ref_dir, out_dir)

# process extracted variants
def process_vars(out_dir):
    for file in os.listdir(out_dir):
        if file.endswith('.filtered.vcf.gz'):
            strip_vcf(out_dir, file)
        else:
            print("No filtered.vcf.gz files found.")

# convert from vcf to plink
def convert_vars(out_dir):
    for file in os.listdir(out_dir):
        if file.endswith('filtered_verified.vcf.gz'):
            vcf_to_bed(out_dir, file)
        else:
            print("No filtered_verfified.vcf files found.")

# run processing functions
def preprocess_vcf(metadata_file, input_dir, ref_dir, out_dir):
    extract_vars(metadata_file, input_dir, ref_dir, out_dir)
    process_vars(out_dir)
    convert_vars(out_dir)

#########################
### merge bed files  ####
#########################

# define function to copy marker file to current wd
def copy_marker(ref_dir, ped_base, output_dir):
    '''
    copy 1000g bed/bim/fam marker files to same dir as vcf files to merge
    :param marker_no_extension: name of marker file without file extension
    :param current_path: vcf_dir to copy marker file to
    :return: copies marker bed/bim/fam file to directory where vcf files are located
    '''
    print "Copying markfer file to {}".format(output_dir)
    copy_marker_cmd = "cp {}/{}.* {}".format(ref_dir, ped_base, output_dir)
    subprocess_cmd(copy_marker_cmd, output_dir)

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

def create_merge_text(input_vcf, merge_list, vcf_dir):
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

def merge_vcf(vcf_dir, ped_base):
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
            create_merge_text(file, merged_file, vcf_dir)
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

def make_pop(ref_dir, fam_dir, merged_fam, kval):
    '''
    Create .pop file with 5 known major populations from 1000g and '-' for unknown from vcf files
    Assumes ref_dir contains tsv files in format subject_id | pop named super_pop.txt and sub_pop.txt
    :param fam_dir path to directory containing merged fam files
    :param ref_dir: path to location of reference files used to create marker set
    :param merged_fam: path to merged fam file to create pop file for
    :param kval: specify 5 for super population and 26 for subpopulations
    :return .pop file for merged bed/bim/fam set as input for supervised admixture
    '''
    # read in _merged.fam file
    in_fam = "{}/{}".format(fam_dir, merged_fam)
    merged_df = pd.read_csv(in_fam, sep=' ', usecols=[1], names=['sample'])
    fam_basename = str('.'.join(merged_fam.split('.')[:-1]) if '.' in merged_fam else merged_fam)
    if kval == 5:
        # pull population from the 1000g map file
        in_file = "{}/super_pop.txt".format(ref_dir)
        ref_map = pd.read_csv(in_file, sep='\t', header=0)
        out_df = pd.merge(merged_df, ref_map, how='left', on = 'sample')
        out_df = out_df['super_pop']
        out_file = "{}/{}.pop".format(fam_dir, fam_basename)
        print ("Saving pop file as {}".format(out_file))
        out_df.to_csv(out_file, sep='\t', na_rep='-', header=False, index=False)
    elif kval == 26:
        # pull population from the 1000g map file
        in_file = "{}/sub_pop.txt".format(ref_dir)
        ref_map = pd.read_csv(in_file, sep='\t', header=0)
        out_df = pd.merge(merged_df, ref_map, how='left', left_on= 'Individual ID', right_on= 'sample')
        out_df = out_df['pop']
        out_file = "{}/{}.pop".format(fam_dir, fam_basename)
        print ("Saving pop file as {}".format(out_file))
        out_df.to_csv(out_file, sep='\t', na_rep='-', header=False, index=False)

def run_admix(admix_path, num_cores, merged_bed, bed_dir, kval):
    '''
    run admixture on input merged bed files
    :param admix_path: path to admixture program
    :param num_cores: number of cores to use for running admixture
    :param merged_bed: bed file of merged marker plus unknown
    :param bed_dir: path to bed file directory
    :param kval: choose either 5 or 26 to run admix with super or sub populations
    :return: P and Q files with admixture results
    '''
    admix_cmd = "nohup {} -j{} {} --supervised {}".format(admixture_path, int(num_cores), merged_bed, kval)
    subprocess_cmd(admix_cmd, bed_dir)

#######################
### Process Results ###
#######################





########################################################

if __name__ == '__main__':

    ### File paths ###

    # path to input vcf files
    vcf_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2'
    # path to reference ped/map files
    ref_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf'
    # ped/map file basename
    ped_base = 'ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
    # path to illumina metadata file
    ill_metadata = "{}/ill_metadata.txt".format(ref_dir)
    # choose either 5 to run admix with super pops, 26 for sub-populations, or 'both'
    pop_kval = 5
    # number of cores to run admixture
    admix_cores = 20

    ### Setup ###

    # make directory to store output
    filtered_out = "{}/filtered/".format(vcf_dir)
    if not os.path.exists(filtered_out):
        os.makedirs(filtered_out)

    # mark true if reference marker needs to be created, else False
    create_marker = 'False'

    ### Run ADMIXTURE ###

    # if create_marker == 'True':
    #     print ("Making marker... ")
    #     # TODO change ref_panel to be determined programatically
    #     make_marker(ref_panel, plink_path, ref_dir, ped_base)

    #preprocess_vcf(ill_metadata, vcf_dir, ref_dir, filtered_out)

    merge_vcf(filtered_out, ped_base)

    # for file in os.listdir(filtered_out):
    #     # check for merge errors and resolve
    #     check_merge_errors(vcf_dir, file)
    #     if file.endswith('_merged.fam'):
    #         # create pop file for each merged fam file
    #         make_pop(ref_dir, filtered_out, file, pop_kval)
    #
    # for file in os.listdir(filtered_out):
    #     if file.endswith('_merged.bed'):
    #         print ("Running admix for {}".format(file))
    #         run_admix(admixture_path, admix_cores, file, filtered_out, pop_kval)




    # read admixture results in
    # q_file = pd.read_csv("/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2/filtered/2016-03-24_102-00144-01.LP6005638-DNA_D06_merged.5.Q", header=None, sep=' ')
    # # read in pop file
    # pop_file = pd.read_csv("/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2/filtered/2016-03-24_102-00144-01.LP6005638-DNA_D06_merged.pop", header=None, sep=' ')
    # # read in fam file
    # fam_file = pd.read_csv("/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2/filtered/2016-03-24_102-00144-01.LP6005638-DNA_D06_merged.fam", header=None, sep= ' ', usecols=[1])
    #
    # print fam_file.head()
    # frames = [fam_file, pop_file, q_file]
    # test = pd.concat(frames, ignore_index=True, keys=None, axis=1)
    # print pd.unique(test[1])
    # print test[(test[1] == 'EAS')].head()
    # print test[(test[1] == '-')]


# TODO add program path as arg to each function
# TODO consistent variable names as input args to functions
# TODO global var with vcf base name and pass to functions
# TODO make basename function and replace this in functions
# TODO all required files in ref_dir and automatically look for file name in func



# TODO each pop file has same subject id, figure out why: error happens when creating metadata.txt