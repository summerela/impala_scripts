##############################################################
## update the following variables before running the script ##
##############################################################
# impala connection strings
impala_host = 'glados18'
impala_port = '21050'

# specify input variants db and table
input_db = 'p7_product'
input_table = 'vars_to_snpeff'

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

####################
## import modules ##
####################
import pandas as pd
from impala.dbapi import connect
import time
import csv
import subprocess

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

## create connection to impala
conn=connect(host=impala_host, port=impala_port, timeout=10000)
cur = conn.cursor()

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

### download variants by row and chromosome
def create_vcf(db_name, table_name, chrom_name):
    # create named file for each chromosome
    vcf_out = 'chr' + chrom_name + '_' + out_name + '.vcf'
    # create header for each chromosome file
    # TODO: make this one function instead of calling header function
    create_header(vcf_out)
    # connect to vars_to_snpeff table
    get_vars = "SELECT chrom, pos, id, ref, alt, qual, filter, info, form, sample from {}.{} WHERE chrom = '{}' order by pos limit 5".format(input_db, input_table, chrom_name)
    cur.execute(get_vars)
    # write variants to file row by row to save memory
    with open(vcf_out, 'a') as csvfile:
        for row in cur:
            writer = csv.writer(csvfile, delimiter="\t", lineterminator = '\n')
            writer.writerow(row)

# # #chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M', 'MT']
chroms = ['1','2','3']

# download each chromosome in input_table and turn into vcf file
for chrom in chroms:
    print "Creating VCF files for chromosome {}... \n".format(chrom)
    create_vcf(input_db, input_table, chrom)

############################################################
# annotate variants with coding consequences using snpeff ##
############################################################
for chrom in chroms:
    print "Annotating coding consequences for chromosome {} with snpeff... \n".format(chrom)
    # create names for input and output files
    vcf_in = 'chr' + chrom + '_' + out_name + '.vcf'
    vcf_out = 'chr' + chrom + '_' + out_name + '_snpeff.vcf'
    # create the file and run snpeff
    with open(vcf_out, "w") as f:
        try:
            subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "-t", "-v", "-noStats", "GRCh37.74", vcf_in], stdout=f)
        except subprocess.CalledProcessError as e:
             print e.output

##########################################################
## Output SnpEff effects as tsv file, one effect per line ##
############################################################
for chrom in chroms:
    print "Parsing snpeff output for chromosome {}... \n".format(chrom)
    vcf_in = 'chr' + chrom + '_' + out_name + '_snpeff.vcf'
    tsv_out = 'chr' + chrom + '_' + out_name + '.tsv'
    # create command to parse snpeff
    snpout_cmd = 'cat {} | {} | {} -jar {} extractFields \
    - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
    "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
    "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(vcf_in, snpeff_oneperline_perl, \
    java_path, snpsift_jar,tsv_out)
    # call subprocess and communicate to pipe output between commands
    ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    print ps.communicate()[0]

####################
## Remove Header  ##
####################
# remove header need for running snpeff to create out own column names on impala
for chrom in chroms:
    print "Removing header for chromosome {} upload to impala... \n".format(chrom)
    tsv_in = 'chr' + chrom + '_' + out_name + '.tsv'
    tsv_out = 'chr' + chrom + '_' + out_name + '_final.tsv'
    tsv_cmd = "sed '1d' {} > {}".format(tsv_in,tsv_out)
    tsv_proc = subprocess.Popen(tsv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print tsv_proc.communicate()[0]

###############################
# ## Upload results to impala ##
# ###############################
import datetime
now = datetime.datetime.now()

# define output path on hdfs
out_path = "{}snpeff_{}".format(hdfs_path, str(now.strftime("%Y%m%d")))

 # put each file in the snpeff directory
for chrom in chroms:
    print "Uploading chromosome {} to HDFS... \n".format(chrom)
    tsv_out = './chr' + chrom + '_' + out_name + '_final.tsv'
    hdfs_cmd = 'hdfs dfs -put {} {}'.format(tsv_out, out_path)
    hdfs_proc = subprocess.Popen(hdfs_cmd, shell=True, stderr=subprocess.STDOUT)
    print hdfs_proc.communicate()[0]

####################################
## Create table to store results  ##
####################################
# drop the table if it already exists
drop_coding = "drop table if exists p7_product.coding_consequences"
cur.execute(drop_coding)

# create empty table to store results
create_coding= '''
create table p7_product.coding_consequences
     (chrom string,
      pos int,
      ref string,
      alt string,
      gene string,
      gene_id string,
      effect string,
      impact string,
      feature string,
      feature_id string,
      biotype string,
      rank int,
      hgvs_c string,
      hgvs_p string)
  row format delimited
  fields terminated by '\t'
'''
cur.execute(create_coding)

##############################
# Insert results into table ##
##############################
# load hdfs files into table
load_query = '''
load data inpath '{}' into table p7_product.coding_consequences
'''.format(out_path)
cur.execute(load_query)

############################
# compute stats on table ##
############################
cur.execute("compute stats  p7_product.coding_consequences")

# TODO add script to make table partitioned

#####################
# close connection ##
#####################
cur.close()



