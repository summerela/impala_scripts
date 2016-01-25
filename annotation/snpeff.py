##############################################################
## update the following variables before running the script ##
##############################################################
# impala connection strings
impala_host = 'glados18'
impala_port = '21050'
impala_user_name = 'selasady'

# specify input variants db and table
input_db = 'p7_product'
input_table = 'dbnsfp_vars'

# prefix for output files
out_name = 'global_vars'

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
chroms = map( str, range(9,23) ) + ['X','Y','M']

##########################################
## create vcf files for each chromosome ##
##########################################
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
    gene_vars = "SELECT chrom, pos, \
    CASE when rs_id is null then '.' \
     else rs_id \
    END AS rs_id,ref, alt, '.' as qual, '.' as filter, '.' as info, '.' as form, '.' as sample \
            from {}.{} WHERE chrom = '{}' order by pos".format(db_name, table_name, chrom_name)
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
    intergenic_vars = "SELECT chrom, pos, \
    CASE when rs_id is null then '.' \
     else rs_id \
    END AS rs_id,ref, alt, '.' as qual, '.' as filter, '.' as info, '.' as form, '.' as sample \
            from {}.{} WHERE chrom = '{}' and gene_name is null order by pos".format(db_name, table_name, chrom_name)
    cur.execute(intergenic_vars)
    int_vars = as_pandas(cur)
    # write variants to file
    if len(int_vars) > 0:
        # create header for each chromosome file
        create_header(vcf_out)
        print "Creating VCF files for chromosome {} intergenic variants... \n".format(chrom_name)
        int_vars.to_csv(vcf_out, sep='\t', index=None, mode='a', header=False)
    else:
        print "No intergenic variants found for chromosome {} \n".format(chrom_name)

# # download each chromosome in input_table and turn into vcf file
for chrom in chroms:
    create_vcf(input_db, input_table, chrom)


# TODO modify pipeline to enable intergenic annotation when snpeff is fixed
# for chrom in chroms:
#     intergenic_vcf(input_db, input_table, chrom)

# ##################################################################
# # check vcf formatting with vcfBareBones.pl from snpeff scripts ##
# ##################################################################
# process all vcf files created from the query
for file in os.listdir(os.getcwd()):
    if any(file.endswith(x) for x in ((out_name + '.vcf'), (out_name + '_intergenic.vcf'))):
        print "Verifying VCF format for {}... \n".format(file)
        vcf_checked_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_verified.vcf'
        # create the file and run snpeff
        with open(vcf_checked_out, "w") as out_file:
            try:
                subprocess.call(['perl', vcf_basic, file], stdout=out_file)
            except subprocess.CalledProcessError as e:
                 print e.output


############################################################
# annotate variants with coding consequences using snpeff ##
############################################################
for file in os.listdir(os.getcwd()):
    # run intergenic variants through snpeff using 'closest' feature to annotate to nearest gene
    if file.endswith('intergenic_verified.vcf'):
        print "Annotating coding consequences for {} with snpeff... \n".format(file)
        # create names for input and output files
        vcf_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_snpeff.vcf'
        # create the file and run snpeff
        with open(vcf_out, "w") as f:
            try:
                subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "closest", "-t", "-v", "GRCh37.75", file], stdout=f)
            except subprocess.CalledProcessError as e:
                 print e.output
    # run non-intergenic variants through snpeff
    elif file.endswith(out_name + '_verified.vcf'):
        print "Annotating coding consequences for {} with snpeff... \n".format(file)
        # create names for input and output files
        vcf_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_snpeff.vcf'
        # create the file and run snpeff
        with open(vcf_out, "w") as f:
            try:
                subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "-t", "-v", "GRCh37.75", file], stdout=f)
            except subprocess.CalledProcessError as e:
                 print e.output

##########################################################
## Output SnpEff effects as tsv file, one effect per line ##
############################################################
for file in os.listdir(os.getcwd()):
    if file.endswith('_snpeff.vcf'):
        print "Parsing snpeff output for {}... \n".format(file)
        tsv_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_parsed.tsv'
        # create command to parse snpeff
        snpout_cmd = 'cat {} | {} | {} -jar {} extractFields \
        - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
        "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].DISTANCE" \
        "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(file, snpeff_oneperline_perl, \
        java_path, snpsift_jar,tsv_out)
        # call subprocess and communicate to pipe output between commands
        ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        print ps.communicate()[0]


# TODO enable function when snpeff closest function is fixed by Broad Inst.

import copy

for file in os.listdir(os.getcwd()):
    if file.endswith('intergenic_verified_snpeff.vcf'):
        print "Parsing snpeff output for {}... \n".format(file)
        # create command to parse snpeff
        df = pd.read_csv(file, sep='\t', skiprows=3, usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO'])

        for row in df.itertuples():
            fixed_cols = list(row[:-1])
            parts = row[-1].split('|')

            # the first one is CLOSEST=0 or something
            fixed_cols.append(parts.pop(0))

            for el in parts:
                my_row = copy.deepcopy(fixed_cols)
                my_row.extend(el.split(','))

                out = ','.join(map(str, my_row))

                with open('./test.csv','ab') as outfile:
                    writer = csv.writer(outfile, lineterminator='\n', newline='')
                    writer.writerow([out])


############################################
## Remove Header and add pos_block column ##
############################################
# remove header that was needed for running snpeff
for file in os.listdir(os.getcwd()):
    if file.endswith('_parsed.tsv'):
        final_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '_final.csv'
        final_df = pd.read_csv(file, sep='\t', skiprows=1, header=None)
        # cant use seq in pandas df slicing
        final_df = final_df[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0]]
        final_df['pos_block'] = final_df[1].div(1000000).astype(int)
        final_df.to_csv(final_out, sep=',', header=False, index=False)

###############################
## Upload results to hdfs  ##
#############################
# TODO add distance column to non-intergenic variants table when fixed
#
import datetime
now = datetime.datetime.now()
today = str(now.strftime("%Y%m%d"))
#
# # define output path on hdfs
out_path = "{}snpeff_{}".format(hdfs_path, today)
mkdir_cmd = "hdfs dfs -mkdir {}".format(out_path)
mkdir_proc = subprocess.Popen(mkdir_cmd, shell=True, stderr=subprocess.STDOUT)
if mkdir_proc.communicate()[0]:
    print "Errors creating HDFS directory: " + mkdir_proc.communicate()[0]

# put each file in the snpeff directory
for file in os.listdir(os.getcwd()):
    if file.endswith('_final.csv'):
        print "Uploading files to HDFS... \n"
        hdfs_cmd = 'hdfs dfs -put {} {}'.format(file, out_path)
        hdfs_proc = subprocess.Popen(hdfs_cmd, shell=True, stderr=subprocess.STDOUT)
        if hdfs_proc.communicate()[0]:
            print "Errors uploading files to HDFS: " + hdfs_proc.communicate()[0]

# set read/write permissions on directory
chown_dir_cmd = "hdfs dfs -chown -R impala:supergroup {}".format(hdfs_path)
chown_proc = subprocess.Popen(chown_dir_cmd, shell=True, stderr=subprocess.STDOUT)
if chown_proc.communicate()[0]:
    print "Errors setting read/write permissions on HDFS directory: " + chown_proc.communicate()[0]

##############################
# Insert results into table ##
##############################
# drop the table if it already exists
drop_coding = "drop table if exists {}.coding_{}".format(input_db, today)
cur.execute(drop_coding)

# create partitioned table
# create_coding_table = '''
# create table {}.coding_{}
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
#       distance int,
#       hgvs_c string,
#       hgvs_p string
#       )
# PARTITIONED BY (chrom string, pos_block int)
#    '''.format(input_db, today)
# cur.execute(create_coding_table)

# load hdfs files into table
# load_query = '''
# load data inpath '{}' into table {}.coding_{}
# '''.format(out_path, input_db, today)
# cur.execute(load_query)

# ############################
# # add pos_block partition ##
# ############################
#
#
# ############################
# # compute stats on table ##
# ############################
# coding_compstats = "compute stats {}.coding_{}".format(input_db, str(now.strftime("%Y%m%d")))
# cur.execute(coding_compstats)
#
# #####################
# # close connection ##
# #####################
# cur.close()
