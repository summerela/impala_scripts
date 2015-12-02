#!/usr/bin/env python

## snpeff_script.py ##
## pipeline for annotating variants with coding consequences using snpeff ##
## input = fill out variables in the first section
## output = impala table of variants with snpeff coding consequence annotations

##############################################
## update the following variables, then run ##
##############################################
# setup impala and hdfs connections
impala_host = 'glados14'
impala_port_number = '21050'
hdfs_host = 'glados20'
hdfs_port_number = '50070'

# specify input variants db and table to annotate with snpeff
input_db = 'p7_product'
input_table = 'all_vars'

# prefix for output files
out_name = 'all_vars'

# path to place file on hdfs
hdfs_path = '/user/selasady/'

# enter full path to java and required jar files
java_path = '/tools/java/jdk1.7/bin/java'
gatk_jar = '/titan/ITMI1/workspaces/users/selasady/tools/GenomeAnalysisTK.jar'
#  enter path to reference fasta used for verifying vcf with gatk
ref_fasta = '/titan/ITMI1/workspaces/users/selasady/tools/human_g1k_v37.fasta'
snpeff_jar = '/titan/ITMI1/workspaces/users/selasady/tools/snpEff/snpEff.jar'
snpeff_oneperline_perl = '/titan/ITMI1/workspaces/users/selasady/tools/snpEff/scripts/vcfEffOnePerLine.pl'
snpsift_jar = '/titan/ITMI1/workspaces/users/selasady/tools/snpEff/SnpSift.jar'

#######################
## create connection ##
#######################
print "Creating a connection to impala.. \n"

import ibis
import os

# connect to impala with ibis
hdfs_port = os.environ.get(hdfs_host, hdfs_port_number)
hdfs = ibis.hdfs_connect(host=hdfs_host, port=hdfs_port, user='hdfs')
con = ibis.impala.connect(host=impala_host, port=impala_port_number, timeout=120)

# enable interactive mode
ibis.options.interactive = True

#################################################
## download required columns and fill in nulls ##
#################################################
import time
import pandas as pd

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

def create_vcf(tbl_name, db_name):
    # create ibis object from distinct vars table
    distinct_vars = con.table(tbl_name, database=db_name) 
    # limit table to just the columns we need to output to vcf
    distinct_df = distinct_vars['chrom', 'pos', 'rs_id', 'ref', 'alt']
    # download table from ibis table connection object to local memory
    distinct_df = distinct_df.execute(limit=100000000000)
    # add blank fields for vcf format
    distinct_df['QUAL'] = '30'
    distinct_df['FILTER'] = 'PASS'
    distinct_df['INFO'] = '.'
    distinct_df['FORMAT'] = 'GT:'
    distinct_df['subject'] = './.'
    # rename chrom column to match vcf header
    distinct_df =distinct_df.rename(columns = {'CHROM':'#CHROM'})
    # rename rs_id column to match vcf header
    distinct_df =distinct_df.rename(columns = {'rs_id':'ID'})
    # uppercase column names to match vcf header
    distinct_df.columns = [x.upper() for x in distinct_df.columns]
    # remove duplicated rows
    distinct_df.drop_duplicates(inplace=True)
    # replace pandas null with '.' for vcf format
    distinct_df['ID'].replace([None], ['.'], inplace=True)
    return distinct_df

# download table from ibis table connection object to local memory
print "Downloading variants... \n"
distinct_df = create_vcf(input_table, input_db)

# display preview of results
print "Variants downloaded... \n" 
print distinct_df.head(5) 
print "\n Total rows = " + str(len(distinct_df))

###################
## format as vcf ##
###################
# create vcf header
def create_header(outfile_name):
   # create vcf header
    lines=[]
    lines.append('##fileformat=VCFv4.0')
    lines.append('##fileDate='+ time.strftime("%y%m%d"))
    lines.append('##reference=grch37 v.74 \n')
    header = '\n'.join(lines)
    out = open(outfile_name, 'wb')
    out.write(header)
    out.close()
    
# create vcf body and append to file with header
def impala_to_vcf(input_df, outfile_name):
    # add blank columns for vcf format and format col names
    input_df.columns = [x.upper() for x in input_df.columns]
    input_df= input_df.rename(columns = {'CHROM':'#CHROM'})    
    # order chromosomes to match ref fastas
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M', 'MT']
    input_df['#CHROM'] = input_df['#CHROM'].astype("category")
    input_df['#CHROM'].cat.set_categories(chroms, inplace=True)
    # sort file by chrom then pos
    input_df = input_df.sort(['#CHROM', 'POS'])
    # write to file for conversion to vcf
    input_df.to_csv(outfile_name, header=True, encoding='utf-8', sep="\t", index=False, mode='a')

vcf_out = out_name + '.vcf'

print "Formatting as VCF file... \n"
create_header(vcf_out)
impala_to_vcf(distinct_df, vcf_out)

####################################################
##  Verify VCF format using GATK ValidateVariants ##
####################################################
import subprocess
print "Validating VCF with GATK ValidateVariants... \n"
subprocess.call([java_path, "-Xmx32g", "-jar", gatk_jar, "-T ", \
 	"ValidateVariants", "-R", ref_fasta, "-V", vcf_out, "--validationTypeToExclude", "ALL"])

##############################################################
##  Annotate variants with coding consequences using snpeff ##
##############################################################
snpeff_out = out_name + '_snpeff.vcf'
print "Running snpeff... \n"
f = open(snpeff_out, "w")
subprocess.call([java_path, "-Xmx32g", "-jar", snpeff_jar, "-t", "-v", "-noStats", "GRCh37.74", vcf_out], stdout=f)

###########################################################
## Output SnpEff effects as tsv file, one effect per line ##
############################################################
tsv_out = out_name + '.tsv'
print "Parsing snpeff output... \n"

# call processes and pipe
snpout_cmd = 'cat {} |{} | {} -jar {} extractFields \
            - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
            "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(snpeff_out, snpeff_oneperline_perl, \
                                                         java_path, snpsift_jar,tsv_out)
ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
output = ps.communicate()[0]
print output
            
######################################
## Remove Header and convert to csv ##
######################################
print "Converting results to csv for upload to impala... \n"
final_out = out_name + '.csv'
csv_cmd = "sed '1d' {} | tr '/\t' ',' > {}".format(tsv_out, final_out)
csv_proc = subprocess.Popen(csv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
csv_proc.communicate()[0]

############################
## upload results to hdfs ##
############################
print "Uploading results to HDFS... \n"
# create directory
if hdfs.exists(hdfs_path):
        hdfs.rmdir(hdfs_path)
        hdfs.mkdir(hdfs_path)
#        hdfs.put(hdfs_path, final_out, verbose=True)

# upload file to hdfs
hdfs.put(hdfs_path, final_out, overwrite=True, verbose=True)

## create table manually until delimited file connector working properly in ibis

############################
## convert csv into table ##
############################
# define table schema for tsv file
#schema = ibis.schema([
#     ('chrom', 'string'), 
#     ('pos', 'int32'),
#     ('ref', 'string'),
#     ('alt', 'string'),
#     ('gene', 'string'),
#     ('gene_id', 'string'),
#     ('effect', 'string'),
#     ('impact', 'string'),
#     ('feature', 'string'),
#     ('feature_id', 'string'),
#     ('biotype', 'string'),
#     ('rank', 'int32'),
#     ('hgvs_c', 'string'),
#     ('hgvs_p', 'string')
#])

#table_name = out_name + '_snpeff'

# create ibis object from  tsv
#print "Creating impala table... \n"
#con.create_table(table_name, con.delimited_file(hdfs_path, schema, delimiter=','), database=input_db, force=True)

print "Snpeff pipeline complete."
