##############################################################
## update the following variables before running the script ##
##############################################################
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

# work mac
java_path = 'java'
gatk_jar =  '/Users/selasady/tools/GenomeAnalysisTK.jar'
ref_fasta = '/Users/selasady/tools/human_g1k_v37.fasta'
snpeff_jar = '/Users/selasady/tools/snpEff/snpEff.jar'
snpeff_oneperline_perl = '/Users/selasady/tools/snpEff/scripts/vcfEffOnePerLine.pl'
snpsift_jar = '/Users/selasady/tools/snpEff//SnpSift.jar'
chrom_splitter = '/Users/selasady/tools/snpEff/scripts/splitChr.pl'

# home pc
# java_path = 'java'
# gatk_jar =  'D:/Documents/tools/GenomeAnalysisTK.jar'
# ref_fasta = 'D:/Documents/tools/human_g1k_v37.fasta'
# snpeff_jar = 'D:/Documents/tools/snpEff/snpEff.jar'
# snpeff_oneperline_perl = 'D;/Documents/tools/snpEff/scripts/vcfEffOnePerLine.pl'
# snpsift_jar = 'D:/Documents/tools/snpEff//SnpSift.jar'
# chrom_splitter = 'D:/Documents/tools/snpEff/scripts/splitChr.pl'

#################################
## create connection to impala ##
#################################
print "Creating a connection to impala.. \n"

import ibis
import os
import sys
import pandas as pd
from impala.dbapi import connect
import time
import csv
import subprocess


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

#################################################
## create vcf files by row for each chromosome ##
#################################################
# # create vcf header
# def create_header(outfile_name):
#    # create vcf header
#     lines=[]
#     lines.append('##fileformat=VCFv4.0')
#     lines.append('##fileDate='+ time.strftime("%y%m%d"))
#     lines.append('##reference=grch37 v.74')
#     lines.append('#CHROM\t' + 'POS\t' + 'ID\t' + 'REF\t' + 'ALT\t' + 'QUAL\t'+ 'FILTER\t' + 'INFO\t' + 'FORMAT\t' + 'SAMPLE\t' + '\n')
#     header = '\n'.join(lines)
#     out = open(outfile_name, 'wb')
#     out.write(header)
#     out.close()
#
# ### download variants by row and chromosome
# def create_vcf(db_name, table_name, chrom_name):
#      vcf_out = 'chr' + chrom_name + '_' + out_name + '.vcf'
#      create_header(vcf_out)
#      # connect to vars_to_snpeff table
#      get_vars = "SELECT chrom, pos, id, ref, alt, qual, filter, info, form, sample from {}.{} WHERE chrom = '{}' order by pos limit 5".format(input_db, input_table, chrom_name)
#      cur.execute(get_vars)
#      for row in cur:
#          with open(vcf_out, 'a') as csvfile:
#              writer = csv.writer(csvfile, delimiter="\t", lineterminator = '\n')
#              writer.writerow(row)
#
# #chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M', 'MT']
chroms = ['1','2','3']
#
# for chrom in chroms:
#     print "Creating VCF files for chromosome {}... \n".format(chrom)
#     create_vcf(input_db, input_table, chrom)
#
# ############################################################
# # annotate variants with coding consequences using snpeff ##
# ############################################################
# for chrom in chroms:
#     print "Annotating coding consequences for chromosome {} with snpeff... \n".format(chrom)
#     vcf_in = 'chr' + chrom + '_' + out_name + '.vcf'
#     vcf_out = 'chr' + chrom + '_' + out_name + '_snpeff.vcf'
#     f = open(vcf_out, "w")
#     try:
#         subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "-t", "-v", "-noStats", "GRCh37.74", vcf_in], stdout=f)
#     except subprocess.CalledProcessError as e:
#         print e.output
#
# ##########################################################
# ## Output SnpEff effects as tsv file, one effect per line ##
# ############################################################
# for chrom in chroms:
#     print "Parsing snpeff output for chromosome {}... \n".format(chrom)
#     vcf_in = 'chr' + chrom + '_' + out_name + '_snpeff.vcf'
#     tsv_out = 'chr' + chrom + '_' + out_name + '.tsv'
#     # call processes and pipe
#     snpout_cmd = 'cat {} | {} | {} -jar {} extractFields \
#     - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
#     "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
#     "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(vcf_in, snpeff_oneperline_perl, \
#     java_path, snpsift_jar,tsv_out)
#     ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
#     ps.communicate()[0]
#
# ####################
# ## Remove Header  ##
# ####################
# for chrom in chroms:
#     print "Removing header for upload to impala... \n"
#     tsv_in = 'chr' + chrom + '_' + out_name + '.tsv'
#     tsv_out = 'chr' + chrom + '_' + out_name + '_final.tsv'
#     tsv_cmd = "sed '1d' {} | tr '/\t' ',' > {}".format(tsv_in,tsv_out)
#     csv_proc = subprocess.Popen(tsv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#     csv_proc.communicate()[0]

# ###############################
# ## Upload results to impala ##
# ###############################
import datetime
now = datetime.datetime.now()

# define output path on hdfs
out_path = "{}/snpeff_{}".format(hdfs_path, str(now.strftime("%Y%m%d")))

# make directory to store output
mkdir_cmd = "hdfs dfs -mkdir {}/snpeff_{}".format(hdfs_path, out_path)
mkdir_proc = subprocess.Popen(mkdir_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
mkdir_proc.communicate()[0]

# put each file in the snpeff directory
for chrom in chroms:
    print "Uploading chromosome {} to HDFS... \n".format(chrom)
    tsv_out = 'chr' + chrom + '_' + out_name + '_final.tsv'
    hdfs_cmd = 'hdfs dfs -put {}/{}'.format(out_path, tsv_out)
    hdfs_proc = subprocess.Popen(hdfs_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    hdfs_proc.communicate()[0]

####################################
## Create table to store results  ##
####################################
# create_coding= '''
# create table p7_product.coding_consequences
#      (pos int,
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
# partitioned by (chrom string)
# '''
# cur.execute(create_coding)

###############################
## Insert results into table ##
###############################
load_query = '''
load data inpath {} into table p7_product.coding_consequences
partition (chrom string)
'''.format(out_path, )
cur.execute(load_query)

#######################
## Close Connections ##
#######################
conn.close()
cur.close()


