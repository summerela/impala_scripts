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
input_table = 'var_test'

# prefix for output files
out_name = 'test_vars'

# path to place file on hdfs
hdfs_path = '/user/selasady/'

java_path = 'java'
gatk_jar =  '/Users/selasady/tools/GenomeAnalysisTK.jar'
ref_fasta = '/Users/selasady/tools/human_g1k_v37.fasta'
snpeff_jar = '/Users/selasady/tools/snpEff/snpEff.jar'
snpeff_oneperline_perl = '/Users/selasady/tools/snpEff/scripts/vcfEffOnePerLine.pl'
snpsift_jar = '/Users/selasady/tools/snpEff//SnpSift.jar'
chrom_splitter = '/Users/selasady/tools/snpEff/scripts/splitChr.pl'

#######################
## create connection ##
#######################
print "Creating a connection to impala.. \n"

import ibis
import os
import sys
from impala.util import as_pandas
import pandas as pd
from impala.dbapi import connect

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

# connect to impala with impyla
conn=connect(host='glados19', port=21050)

# connect to impala with ibis
hdfs_port = os.environ.get(hdfs_host, hdfs_port_number)
hdfs = ibis.hdfs_connect(host=hdfs_host, port=hdfs_port, user='hdfs')
con = ibis.impala.connect(host=impala_host, port=impala_port_number, timeout=120)

# enable interactive mode
ibis.options.interactive = True

### download variants by chromosome
def create_vcf(tbl_name, db_name, chrom_name):
    # create ibis object from distinct vars table
    distinct_vars = con.table(tbl_name, database=db_name)
    # limit table to just the columns we need to output to vcf
    distinct_df = distinct_vars['chrom', 'pos', 'rs_id', 'ref', 'alt']
    # limit table by chromosome
    distinct_df = distinct_df[distinct_df.chrom == chrom_name].limit(5)
    # download table from ibis table connection object to local memory
    for row in distinct_df.execute().iteritems():
        # add blank fields for vcf format
        distinct_df['QUAL'] = '.'
        distinct_df['FILTER'] = '.'
        distinct_df['INFO'] = '.'
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

###################
## format as vcf ##
###################
# import time
#
# # create vcf header
# def create_header(outfile_name):
#    # create vcf header
#     lines=[]
#     lines.append('##fileformat=VCFv4.0')
#     lines.append('##fileDate='+ time.strftime("%y%m%d"))
#     lines.append('##reference=grch37 v.74 \n')
#     header = '\n'.join(lines)
#     out = open(outfile_name, 'wb')
#     out.write(header)
#     out.close()
#
# # create vcf body and append to file with header
# def impala_to_vcf(input_df, outfile_name):
#     # add blank columns for vcf format and format col names
#     input_df.columns = [x.upper() for x in input_df.columns]
#     input_df= input_df.rename(columns = {'CHROM':'#CHROM'})
#     # write to file for conversion to vcf
#     input_df.to_csv(outfile_name, header=True, encoding='utf-8', sep="\t", index=False, mode='a')
#

chroms = ['1','2','3']
for chrom in chroms:
    print "Downloading variants in chromosome {}".format(chrom)
    distinct_df = create_vcf(input_table, input_db, chrom)
#     vcf_out = 'chr' + chrom + '_' + out_name + '.vcf'
#     print "Formatting chromsome {} as VCF file... \n".format(chrom)
#     create_header(vcf_out)
#     impala_to_vcf(distinct_df, vcf_out)
#     print "Finished with test"
    # import subprocess
    # print "Validating VCF with GATK ValidateVariants... \n"
    # subprocess.call([java_path, "-Xmx32g", "-jar", gatk_jar, "-T ", \
 	# "ValidateVariants", "-R", ref_fasta, "-V", vcf_out, "--validationTypeToExclude", "ALL"])
    # snpeff_out =  'chr' + chrom + '_snpeff.vcf'
    # print "Running snpeff... \n"
    # f = open(snpeff_out, "w")
    # subprocess.call([java_path, "-Xmx32g", "-jar", snpeff_jar, "-t", "-v", "-noStats", "GRCh37.74", vcf_out], stdout=f)




