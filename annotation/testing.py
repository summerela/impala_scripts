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
input_table = 'vars_to_snpeff'

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
import pandas as pd
from impala.dbapi import connect
import time
import csv

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

# connect to impala with impyla
conn=connect(host='glados19', port=21050)
cur = conn.cursor()

# connect to impala with ibis
hdfs_port = os.environ.get(hdfs_host, hdfs_port_number)
hdfs = ibis.hdfs_connect(host=hdfs_host, port=hdfs_port, user='hdfs')
con = ibis.impala.connect(host=impala_host, port=impala_port_number, timeout=120)

# enable interactive mode
ibis.options.interactive = True

# create vcf header
def create_header(outfile_name):
   # create vcf header
    lines=[]
    lines.append('##fileformat=VCFv4.0')
    lines.append('##fileDate='+ time.strftime("%y%m%d"))
    lines.append('##reference=grch37 v.74 \n')
    lines.append('CHROM \t' + 'POS\t' + 'ID\t' + 'REF\t' + 'ALT\t' + 'QUAL\t'+ 'FILTER\t' + 'INFO\t' + '\n')
    header = '\n'.join(lines)
    out = open(outfile_name, 'wb')
    out.write(header)
    out.close()

### download variants by row and chromosome
def create_vcf(db_name, table_name, chrom_name):
     vcf_out = 'chr' + chrom_name + '_' + out_name + '.vcf'
     create_header(vcf_out)
     # connect to vars_to_snpeff table
     get_vars = "SELECT chrom, pos, id, ref, alt, qual, filter, info from {}.{} WHERE chrom = '{}'".format(input_db, input_table, chrom_name)
     cur.execute(get_vars)
     for row in cur:
         with open(vcf_out, 'a') as csvfile:
             writer = csv.writer(csvfile, delimiter="\t")
             writer.writerow(row)

chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M', 'MT']
for chrom in chroms:
    create_vcf(input_db, input_table, chrom)






# chroms = ['1','2','3']
# for chrom in chroms:
#     print "Downloading variants in chromosome {}".format(chrom)
#     distinct_df = create_vcf(input_table, input_db, chrom)
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




