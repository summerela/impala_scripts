'''
Run ADMIXTURE
- create reference marker file
- convert VCF of interest to plink bed format
- run admixture
- output results to impala table
'''

####################
### Script Setup ###
####################

# import modules
import os
import subprocess as sp
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

# make directory to store output
def make_out_dir(vcf_dir):
    filtered_out = "{}/filtered/".format(vcf_dir)
    if not os.path.exists(filtered_out):
        os.makedirs(filtered_out)

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

# TODO convert get_base_name to class and document
def get_base_name(in_file, num_to_remove):
    basename = str('.'.join(in_file.split('.')[:-int(num_to_remove)]) if '.' in in_file else in_file)
    return basename

def get_ped_base(ref_dir):
    for file in os.listdir(ref_dir):
        if file.endswith('.ped'):
            ped_base = get_base_name(file, 1)
            return str(ped_base)

def get_gzVCF_base(vcf_in):
    vcf_base = get_base_name(vcf_in, 2)
    return vcf_base

def get_outdir(vcf_dir):
    out_dir = "{}/filtered/".format(vcf_dir)
    return out_dir

####################################
### Create reference marker file ###
####################################

# create text file with subject id's grouped by population for sorting subjects by population
def create_sorter_and_pop(input_dir):
    '''
    :param input_dir: ref_dir containing ref panel file (three columns, tsv, col headers = sample | pop | super_pop )
    :return sorted_ref.txt a two column tsv ( sample | pop ) to be fed to plink for sorting by pop
    '''
    print("\n Creating sorted panel file and pop files for super and sub populations... \n")
    ref_panel = ''.join([f for f in os.listdir(input_dir) if f.endswith('.panel')])
    ref_panel_file = "{}/{}".format(input_dir, str(ref_panel))
    out_panel = "{}/sorted_panel.txt".format(input_dir)
    print("Saving sorted panel file as {}".format(out_panel))
    panel = pd.read_csv(ref_panel_file, sep='\t', usecols=['sample', 'pop', 'super_pop', 'gender'])
    panel.sort_values(by=['super_pop', 'pop'], axis = 0, inplace=True)
    # create super_pop.txt file for super populations
    super_out = panel[['sample', 'pop']]
    super_out.to_csv("{}/super_pop.txt".format(input_dir), index=False, names=['sample', 'super_pop'])
    # create sub_pop.txt for subpopulations
    sub_out = panel[['sample', 'super_pop']]
    sub_out.to_csv("{}/sub_pop.txt".format(input_dir), index=False, names=['sample', 'pop'])
    # create sorted_ref.txt for sorting marker tsv by pop in format family_id | subject_id
    sorter_out = panel[['sample', 'sample']]
    sorter_out.to_csv("{}/sorted_panel.txt".format(input_dir), index=False, sep='\t', header=False)


def create_marker(input_dir, plink):
    '''
    Runs commands to create marker file if make_marker == True
    :param plink: path to plink1.9
    :param ref_dir: path to dir containing plink ped/map files
    :return: marker file to merge with input VCF files for running ADMIXTURE
    '''
    # get basename of ped file in ref_dir
    print ("test")
    # ped_base = get_ped_base(input_dir)
    # filter reference ped file to retain MAF between .05 and 0.5 and LD 50 5 0.5
    # plink_cmd = "{} --file {}/{} --out {}/{}_maf_ld --maf 0.05 --max-maf .49 --indep-pairwise 50 5 0.5".format(plink, input_dir, ped_base, input_dir, ped_base)
    # # retain only maf/ld pruned variants and convert ped to bed/bim/fam file using plink
    # pruned_file = '{}/{}_maf_ld.prune.in'.format(input_dir, ped_base)
    # ped2bed_cmd = '{} --file {}/{} --extract {} --make-bed --indiv-sort file {}/sorted_ref.txt --out {}/{}'.format(plink, input_dir, ped_base, pruned_file, input_dir, input_dir, ped_base)
    # # create list of snps from marker file to extract from input vcf files
    # ref_snps = "{} --bfile {} --write-snplist --out {}_snps".format(plink, ped_base,ped_base)
    # # run the vcf procesing commands
    # marker_cmds = [plink_cmd,ped2bed_cmd, ref_snps]
    # for cmd in marker_cmds:
    #     subprocess_cmd(cmd, input_dir)

# def process_reference(input_dir, plink_path):
    # create_sorter_and_pop(input_dir)
    # create_marker(input_dir, plink_path)

###############################
### Pre-process VCF file(s) ###
###############################

def trim_vcf(input_dir, input_vcf, vcf_verify_path):
    '''
    process all vcf.gz files to retain only chrom, pos, ref, alt and gt
    :param input_dir: path to directory where vcf files are located
    :param input_vcf: vcf.gz file to trim
    :return: _trimmed.vcf.gz files to merge with marker file before running ADMIXTURE
    '''
    print ("\n Trimming VCF file {}... \n".format(input_vcf))
    vcf_base = get_gzVCF_base(input_vcf)
    vcf_out = "{}_trimmed.vcf.gz".format(vcf_base)
    out_dir = get_outdir(input_dir)
    snp_verify_cmd = 'nohup zcat {}/{} | {} | gzip > {}/{} '.format(input_dir, input_vcf, vcf_verify_path ,out_dir, vcf_out)
    subprocess_cmd(snp_verify_cmd, input_dir)

def vcf_to_bed(input_dir, input_vcf):
    '''
    convert stripped vcf to bed/bim/fam with plink 1.9
    :param input_dir: directory containing _verified.vcf.gz created from strip_vcf()
    :param input_vcf: vcf file to process
    :return: stripped, verified vcf files converted to plink binary (bed) format
    '''
    bed_out = get_gzVCF_base(input_vcf)
    out_dir = get_outdir(input_dir)
    vcf2plink_cmd = "nohup {} --vcf {}/{} --double-id --biallelic-only strict --memory 300000 --geno 0.1 --allow-no-sex --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out {}/{}".format(plink_path, input_dir, input_vcf, out_dir, bed_out)
    print ("Converting {} from vcf to bed/bim/fam").format(input_vcf)
    subprocess_cmd(vcf2plink_cmd, input_dir)

# subset vcf for maker regions while converting to plink binary format
def subset_bed(input_dir, input_vcf, ref_snps):
    '''
    subset vcf file using plink output snps from marker set
    :param input_vcf: vcf file to subset. assumes .vcf.gz
    :param ref_snps: plink output of snp list from marker file using make_marker.ref_snps command
    :return: metadata subset dataframe containing vcf file info to feed to update_meta_path()
    '''
    vcf_basename = get_gzVCF_base(vcf_dir)
    out_dir = get_outdir(input_dir)
    subset_bed_cmd = "nohup {} --bfile {} --extract {} --make-bed --out {}/{}".format(plink_path,input_vcf,ref_snps, out_dir,vcf_basename)
    subprocess_cmd(subset_bed_cmd, input_dir)

def get_ref_snps(ref_dir):
    snp_list = [f for f in os.listdir(ref_dir) if f.endswith('.snplist')]
    return snp_list

# run processing functions
def run_trim(input_dir, ref_dir, vcf_verify_path):
    for file in os.listdir(input_dir):
        if file.endswith('.vcf.gz'):
            trim_vcf(input_dir, file, vcf_verify_path)

def convert_vcf(input_dir, ref_dir):
    if file.endswith('_trimmed.vcf.gz'):
        vcf_to_bed(input_dir, file)

def run_subset(input_dir):
    for file in os.listdir(input_dir):
        if file.endswith('.bed'):
            ref_snps = get_ref_snps(ref_dir)
            subset_bed(input_dir, file, ref_snps)

def process_vcf(input_dir, ref_dir, vcf_verify_path):
    run_trim(input_dir, ref_dir, vcf_verify_path)
    convert_vcf(input_dir, ref_dir)
    run_subset(input_dir)



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
    #subprocess_cmd(merge_cmd, vcf_dir)
    print merge_cmd

def merge_vcf(vcf_dir, ref_dir, ped_base):
    '''
    run pipeline to merge bed files with marker set
    :param vcf_dir: directory path to vcf files to merge
    :return: merged bed/bim/fam file to run through ADMIXTURE program
    '''
    # copy marker file to vcf_dir
    copy_marker(ref_dir, ped_base, vcf_dir)
    for file in os.listdir(vcf_dir):
        if file.endswith('_verified.bed'):
            #merged_file = findBedStems(vcf_dir, file, ped_base)
            #create_merge_text(file, merged_file, vcf_dir)
            plink2_merge(file, vcf_dir)

#######################
###  run admixture  ###
#######################

def make_pop(ref_dir, out_dir, merged_fam, kval):
    '''
    Create .pop file with 5 known major populations from 1000g and '-' for unknown from vcf files
    Assumes ref_dir contains tsv files in format subject_id | pop named super_pop.txt and sub_pop.txt
    :param out_dir path to directory containing merged fam files
    :param ref_dir: path to location of reference files used to create marker set
    :param merged_fam: path to merged fam file to create pop file for
    :param kval: specify 5 for super population and 26 for subpopulations
    :return .pop file for merged bed/bim/fam set as input for supervised admixture
    '''
    # read in _merged.fam file
    in_fam = "{}/{}".format(out_dir, merged_fam)
    merged_df = pd.read_csv(in_fam, sep=' ', usecols=[1], names=['sample'])
    fam_basename = str('.'.join(merged_fam.split('.')[:-1]) if '.' in merged_fam else merged_fam)
    if kval == 5:
        # pull population from the 1000g map file
        in_file = "{}/super_pop.txt".format(ref_dir)
        ref_map = pd.read_csv(in_file, sep='\t', header=0)
        out_df = pd.merge(merged_df, ref_map, how='left', on = 'sample')
        out_df = out_df['super_pop']
        out_file = "{}/{}.pop".format(out_dir, fam_basename)
        print ("Saving pop file as {}".format(out_file))
        out_df.to_csv(out_file, sep='\t', na_rep='-', header=False, index=False)
    elif kval == 26:
        # pull population from the 1000g map file
        in_file = "{}/sub_pop.txt".format(ref_dir)
        ref_map = pd.read_csv(in_file, sep='\t', header=0)
        out_df = pd.merge(merged_df, ref_map, how='left', left_on= 'Individual ID', right_on= 'sample')
        out_df = out_df['pop']
        out_file = "{}/{}.pop".format(out_dir, fam_basename)
        print ("Saving pop file as {}".format(out_file))
        out_df.to_csv(out_file, sep='\t', na_rep='-', header=False, index=False)

def admix_cmd(admix_path, num_cores, merged_bed, out_dir, kval):
    '''
    run admixture on input merged bed files
    :param admix_path: path to admixture program
    :param num_cores: number of cores to use for running admixture
    :param merged_bed: bed file of merged marker plus unknown
    :param out_dir: path to bed file directory
    :param kval: choose either 5 or 26 to run admix with super or sub populations
    :return: P and Q files with admixture results
    '''
    admix_cmd = "nohup {} -j{} {} --supervised {}".format(admixture_path, int(num_cores), merged_bed, kval)
    subprocess_cmd(admix_cmd, out_dir)

def run_pop(ref_dir, out_dir, kval):
    for file in os.listdir(out_dir):
        if file.endswith('_merged.fam'):
            make_pop(ref_dir, out_dir, file, kval)

def run_admix(admix_path, num_cores, out_dir, kval):
    for file in os.listdir(out_dir):
        if file.endswith('_merged.bed'):
            admix_cmd(admix_path, num_cores, file, out_dir, kval)

#######################
### Process Results ###
#######################





########################################################

if __name__ == '__main__':

    ### File paths ###

    # path to input vcf files
    vcf_file_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test2/'
    # path to reference ped/map files
    ref_file_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'
    # ped/map file basename
    ped_base = 'ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
    # choose either 5 to run admix with super pops, 26 for sub-populations, or 'both'
    pop_kval = 5
    # number of cores to run admixture
    admix_cores = 20

    # mark true if reference marker needs to be created, else False
    create_marker = 'True'

    ### Run ADMIXTURE ###

    make_out_dir(vcf_file_dir)
    filtered_out = "{}/filtered/".format(vcf_file_dir)

    if create_marker == 'True':
        print ("Making marker... ")
        #create_sorter_and_pop(ref_file_dir)
        create_marker(ref_file_dir, plink_path)





    #merge_vcf(filtered_out, ped_base)
    #run_pop(ref_dir, filtered_out, pop_kval)
    #run_admix(admixture_path, admix_cores, filtered_out, pop_kval)


    # for file in os.listdir(filtered_out):
    #     if file.endswith(".Q"):
    #         # read admixture results in
    #         q_file = pd.read_csv("{}/{}".format(filtered_out, file), header=None, sep=' ')
    #         base_start = str('.'.join(file.split('.')[:-2]) if '.' in file else file)
    #         pop_file = "{}/{}.pop".format(filtered_out, base_start)
    #         # read in pop file
    #         pop_file = pd.read_csv(pop_file, header=None, sep=' ')
    #         # read in fam file
    #         fam_file = "{}/{}.fam".format(filtered_out, base_start)
    #         fam_file = pd.read_csv(fam_file, header=None, sep= ' ', usecols=[1])
    #         frames = [fam_file, pop_file, q_file]
    #         out_frame = pd.concat(frames, ignore_index=True, keys=None, axis=1)
    #         out_frame.columns= ['subject', 'known_pop', 'EUR', 'EAS', 'AMR', 'SAS', 'AFR']
    #         print out_frame[(out_frame['known_pop'] == '-')]




# TODO add program path as arg to each function
# TODO consistent variable names as input args to functions
# TODO global var with vcf base name and pass to functions
# TODO all required files in ref_dir and automatically look for file name in func


##### unused code snippets ###

# # # convert the binary bed file to a readable map file to get genomic locations for extracting vcf regions
# # map_pos_cmd = "{} --noweb --bfile {}/{} --recode 12 --out {}/{}_map".format(plink_path, ref_dir, ped_base, ref_dir, ped_base)
# # convert bed file to vcf
# bed_to_vcf_cmd = "{} --file {}  --recode vcf --out {}".format(plink_path, ped_base, ped_base)
# # extract chrom pos ref from marker vcf
# make_region = r'''cat {}.vcf  | grep -v '^#' | awk '{{OFS="\t";print "chr"$1, $2, $4}}' | gzip > {}_markers.gz'''.format(ped_base, ped_base)





# def update_meta_path(vcf_subset, input_dir, input_vcf, out_dir):
#     '''
#     Update vcf file path for subset metadata file for current vcf location
#     :param vcf_subset: dataframe subset of platform metadata file created from subset_metadata()
#     :param input_dir: directory containing vcf file
#     :param input_vcf: vcf file to subset metadata for and update the file path to location of this vcf
#     :param out_dir: directory to store output metadata text
#     :return: subset metadata file containing information for input_vcf as tsv
#     '''
#     vcf_path = os.path.abspath(os.path.join(input_dir, input_vcf))
#     print ("\n Creating metadata file for {}".format(vcf_path))
#     vcf_subset['VCF'] = vcf_path
#     vcf_basename = str('.'.join(input_vcf.split('.')[:1]) if '.' in input_vcf else input_vcf)
#     out_file = "{}/{}_metadata.txt".format(out_dir, vcf_basename)
#     vcf_subset.to_csv(out_file, sep='\t', index=False, na_rep="NA")
#
# # extract regions from vcf file that match marker regions
# def extract_variants(input_vcf, ped_base, ref_dir, out_dir):
#     '''
#     Runs extract_regions_stream_job.pl written by Denise Maulden to extract regions from input vcf files
#      that match regions in the marker file to save computational time
#     :param input_vcf: path to vcf file to extract regions from
#     :param ped_base: base name of input ped file used to create marker file
#     :param ref_dir: dictory path to reference marker file
#     :param metadata_file: file with metadata and location of input vcf files made from create_metadata_file()
#     :param out_dir: directory to write output vcf's
#     :return: *_filtered.vcf.gz files retaining only regions that match the maker file regions
#     '''
#     print ("Extracting marker regions from vcf files... \n")
#     vcf_basename = str('.'.join(input_vcf.split('.')[:1]) if '.' in input_vcf else input_vcf)
#     region_file = "{}/{}_markers.gz".format(ref_dir, ped_base)
#     metadata_file = "{}{}_metadata.txt".format(out_dir, vcf_basename)
#     extract_cmd = "/tools/bin/perl {} --regionFile {} --metadata {} --outDir {}  --compressionType bzip2".format(extract_script, region_file, metadata_file, out_dir)
#     subprocess_cmd(extract_cmd, vcf_dir)


#meta_names = ['Vendor',  'Bucket',  'Study',   'Family',  'Genome',  'Subject', 'Sample',  'Assembly', 'Gestalt ID', 'Gender',  'Term'
                  # 'Country of Birth', 'Member',  'VCF', 'PATH', 'S3','BIGDATA', 'GESTALT',  'ITMI_MAPPINGS', 'COUNT', 'S3URL', 'BigdataURL',      'GestaltURL', 'MappingsURL']
    # read in metadata file
# meta_file = pd.read_csv(metadata_file, sep='\t', header=0)
# # subset metadata file and update vcf file paths
# vcf_subset = meta_file[(meta_file['Genome'].str.contains(vcf_basename))]
# return vcf_subset