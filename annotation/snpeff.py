##############################################################
## update the following variables before running the script ##
##############################################################
# impala connection strings
impala_host = 'glados18'
impala_port = '21050'
impala_user_name = 'selasady'

# specify input variants db and table
input_db = 'p7_product'
input_table = 'test_vars'

# prefix for output files
out_name = 'test_vars'

# home path to user's hdfs directory
hdfs_path = '/user/selasady/'

# unix file paths
java_path = '/tools/java/jdk1.7/bin/java'
gatk_jar =  '/users/selasady/my_titan_itmi/tools/GenomeAnalysisTK.jar'
ref_fasta = '/users/selasady/my_titan_itmi/tools/human_g1k_v37.fasta'
snpeff_jar = '/users/selasady/my_titan_itmi/tools/snpEff/snpEff.jar'
snpeff_oneperline_perl = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/vcfEffOnePerLine.pl'
snpsift_jar = '/users/selasady/my_titan_itmi/tools/snpEff//SnpSift.jar'
chrom_splitter = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/splitChr.pl'
vcf_basic = '/users/selasady/my_titan_itmi/impala_scripts/annotation/parse_vcf.pl'

####################
## import modules ##
####################
import pandas as pd
from impala.dbapi import connect
import time
import csv
import subprocess
import numpy as np
from impala.util import as_pandas
import os

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

## create connection to impala
conn=connect(host=impala_host, port=impala_port, timeout=10000, user=impala_user_name)
cur = conn.cursor()

# create list of chromosomes to process
chroms = map( str, range(1,23) ) + ['X','Y','M']

#################################################
## create vcf files by row for each chromosome ##
#################################################
# create vcf header
def create_header(outfile_name):
   # create vcf header
    lines=[]
    lines.append('##fileformat=VCFv4.0')
    lines.append('##fileDate='+ time.strftime("%y%m%d"))
    lines.append('##reference=grch37 v.74')
    lines.append('#CHROM\t' + 'POS\t' + 'ID\t' + 'REF\t' + 'ALT\t' + 'QUAL\t'+ 'FILTER\t' + 'INFO\t' + 'FORMAT\t' + 'SAMPLE\t' + '\n')
    header = '\n'.join(lines)
    out = open(outfile_name, 'wb')
    out.write(header)
    out.close()

### download variants that are not intergenic
def create_vcf(db_name, table_name, chrom_name):
    # create named file for each chromosome
    vcf_out = 'chr' + chrom_name + '_' + out_name + '.vcf'
    # connect to vars_to_snpeff table
    gene_vars = "SELECT chrom, pos, rs_id, ref, alt, '.' as qual, '.' as filter, '.' as info, '.' as form, '.' as \
                sample from {}.{} WHERE chrom = '{}' and gene_name is not null order by pos".format(db_name, table_name, chrom_name)
    cur.execute(gene_vars)
    vars = as_pandas(cur)
    # write variants to file
    if len(vars) > 0:
        # create header for each chromosome file
        # TODO: make this one function instead of calling header function
        create_header(vcf_out)
        print "Creating VCF files for chromosome {}... \n".format(chrom_name)
        vars.to_csv(vcf_out, sep='\t', index=None, mode='a', header=False)
    else:
        print "No variants found for chromosome {} \n".format(chrom_name)


### download variants that are in intergenic regions
def intergenic_vcf(db_name, table_name, chrom_name):
    # create named file for each chromosome
    vcf_out = 'chr' + chrom_name + '_' + out_name + '_intergenic.vcf'
    # connect to vars_to_snpeff table
    intergenic_vars = "SELECT chrom, pos, rs_id, ref, alt, '.' as qual, '.' as filter, '.' as info, '.' as form, '.' as \
                sample from {}.{} WHERE chrom = '{}' and gene_name is null order by pos".format(db_name, table_name, chrom_name)
    cur.execute(intergenic_vars)
    int_vars = as_pandas(cur)
    # write variants to file
    if len(int_vars) > 0:
        # create header for each chromosome file
        # TODO: make this one function instead of calling header function
        create_header(vcf_out)
        print "Creating VCF files for chromosome {} intergenic variants... \n".format(chrom_name)
        int_vars.to_csv(vcf_out, sep='\t', index=None, mode='a', header=False)
    else:
        print "No intergenic variants found for chromosome {} \n".format(chrom_name)

# download each chromosome in input_table and turn into vcf file
# for chrom in chroms:
#     create_vcf(input_db, input_table, chrom)
#
# for chrom in chroms:
#     intergenic_vcf(input_db, input_table, chrom)

# ##################################################################
# # check vcf formatting with vcfBareBones.pl from snpeff scripts ##
# ##################################################################
# process all vcf files created from the query
# for file in os.listdir(os.getcwd()):
#     if any(file.endswith(x) for x in ((out_name + '.vcf'), (out_name + '_intergenic.vcf'))):
#         print "Verifying VCF format for {}... \n".format(file)
#         vcf_checked_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_verified.vcf'
#         # create the file and run snpeff
#         with open(vcf_checked_out, "w") as out_file:
#             try:
#                 subprocess.call(['perl', vcf_basic, file], stdout=out_file)
#             except subprocess.CalledProcessError as e:
#                  print e.output

# ############################################################
# # annotate variants with coding consequences using snpeff ##
# ############################################################
# for file in os.listdir(os.getcwd()):
#     # run intergenic variants through snpeff using 'closest' feature to annotate to nearest gene
#     if file.endswith('intergenic_verified.vcf'):
#         print "Annotating coding consequences for {} with snpeff... \n".format(file)
#         # create names for input and output files
#         vcf_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_intergenic_snpeff.vcf'
#         # create the file and run snpeff
#         with open(vcf_out, "w") as f:
#             try:
#                 subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "closest", "-t", "-v", "GRCh37.75", file], stdout=f)
#             except subprocess.CalledProcessError as e:
#                  print e.output
#     # run non-intergenic variants through snpeff
#     elif file.endswith(out_name + '_verified.vcf'):
#         print "Annotating coding consequences for {} with snpeff... \n".format(file)
#         # create names for input and output files
#         vcf_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_snpeff.vcf'
#         # create the file and run snpeff
#         with open(vcf_out, "w") as f:
#             try:
#                 subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "-t", "-v", "GRCh37.75", file], stdout=f)
#             except subprocess.CalledProcessError as e:
#                  print e.output

# ##########################################################
# ## Output SnpEff effects as tsv file, one effect per line ##
# ############################################################
# for file in os.listdir(os.getcwd()):
#     if file.endswith('_snpeff.vcf'):
#         print "Parsing snpeff output for {}... \n".format(file)
#         tsv_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_parsed.tsv'
#         # create command to parse snpeff
#         snpout_cmd = 'cat {} | {} | {} -jar {} extractFields \
#         - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
#         "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
#         "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(file, snpeff_oneperline_perl, \
#         java_path, snpsift_jar,tsv_out)
#         # call subprocess and communicate to pipe output between commands
#         ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
#         print ps.communicate()[0]

for file in os.listdir(os.getcwd()):
    if file.endswith('intergenic_snpeff.vcf'):
        print "Parsing snpeff output for {}... \n".format(file)
        tsv_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_parsed.tsv'
        # create command to parse snpeff
        df = pd.read_csv(file, sep='\t', skiprows=4, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'CLOSEST'])
        print df

#


# ####################
# ## Remove Header  ##
# ####################
# remove header need for running snpeff to create out own column names on impala
# for file in os.listdir(os.getcwd()):
#     if file.endswith('_final.tsv'):
#         print "Removing header for {}... \n".format(file)
#         tsv_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_nohead.tsv'
#         tsv_cmd = "sed '1d' {} > {}".format(file,tsv_out)
#         tsv_proc = subprocess.Popen(tsv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#         print tsv_proc.communicate()[0]
#
# ###############################
# ## Upload results to hdfs  ##
# #############################
# import datetime
# now = datetime.datetime.now()
#
# # define output path on hdfs
# out_path = "{}snpeff_{}".format(hdfs_path, str(now.strftime("%Y%m%d")))
# mkdir_cmd = "hdfs dfs -mkdir {}".format(out_path)
# mkdir_proc = subprocess.Popen(mkdir_cmd, shell=True, stderr=subprocess.STDOUT)
# print mkdir_proc.communicate()[0]
#
#  # put each file in the snpeff directory
# print "Uploading files to HDFS... \n"
# hdfs_cmd = 'hdfs dfs -put chr*_final.tsv {}'.format(out_path)
# hdfs_proc = subprocess.Popen(hdfs_cmd, shell=True, stderr=subprocess.STDOUT)
# print hdfs_proc.communicate()[0]
#
# # set read/write permissions on directory
# chown_dir_cmd = "hdfs dfs -chown -R impala:supergroup {}".format(hdfs_path)
# chown_proc = subprocess.Popen(chown_dir_cmd, shell=True, stderr=subprocess.STDOUT)
# print chown_proc.communicate()[0]
#
# # ####################################
# # ## Create table to store results  ##
# # ####################################
# # drop the table if it already exists
# drop_coding = "drop table if exists {}.coding_consequences".format(input_db)
# cur.execute(drop_coding)
#
# # create empty table to store results
# create_coding= '''
# create table {}.coding_consequences
#      (chrom string,
#       pos int,
#       id string,
#       ref string,
#       alt string,
#       gene string,
#       gene_id string,
#       effect string,
#       impact string,
#       feature string,
#       feature_id string,
#       biotype string,
#       rank int,
#       hgvs_c string,
#       hgvs_p string)
#   row format delimited
#   fields terminated by '\t'
# '''.format(input_db)
# cur.execute(create_coding)
#
# ##############################
# # Insert results into table ##
# ##############################
# # load hdfs files into table
# load_query = '''
# load data inpath '{}' into table {}.coding_consequences
# '''.format(out_path, input_db)
# cur.execute(load_query)
#
# ############################
# # compute stats on table ##
# ############################
# coding_compstats = "compute stats {}.coding_consequences".format(input_db)
# cur.execute(coding_compstats)
#
# ##########################
# # partition final table ##
# ##########################
#
# # drop the table if it already exists
# drop_coding = "drop table if exists {}.all_coding".format(input_db)
# cur.execute(drop_coding)
#
# # create partitioned table
# create_parition = '''
# create table {}.all_coding
#     (
#       pos int,
#       id string,
#       ref string,
#       alt string,
#       gene string,
#       gene_id string,
#       effect string,
#       impact string,
#       feature string,
#       feature_id string,
#       biotype string,
#       rank int,
#       hgvs_c string,
#       hgvs_p string
#       )
# partitioned by (chrom string)'''.format(input_db)
# cur.execute(create_parition)
#
# # insert unpartitioned table
# insert_partition = '''
# insert into table {}.all_coding partition(chrom)
#     (
#       select  pos, id, ref, alt, gene, gene_id, effect, impact, feature, feature_id,
#       biotype, rank, hgvs_c, hgvs_p, chrom
#       from p7_product.coding_consequences
#       )'''.format(input_db)
# cur.execute(insert_partition)
#
# # compute stats
# final_compstats = "compute stats {}.all_coding".format(input_db)
# cur.execute(final_compstats)
#
#
# #####################
# # close connection ##
# #####################
# cur.close()
#
#
#
