# specify table to check for new variants
input_db = "p7_product"
input_table =  "test"
hdfs_path = "/user/selasady/"
result_name = "new_vars"

# system paths to required executables
java_path = '/tools/java/jdk1.7/bin/java'
gatk_jar =  '/users/selasady/my_titan_itmi/tools/GenomeAnalysisTK.jar'
ref_fasta = '/users/selasady/my_titan_itmi/tools/human_g1k_v37.fasta'
snpeff_jar = '/users/selasady/my_titan_itmi/tools/snpEff/snpEff.jar'
snpeff_oneperline_perl = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/vcfEffOnePerLine.pl'
snpsift_jar = '/users/selasady/my_titan_itmi/tools/snpEff//SnpSift.jar'
chrom_splitter = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/splitChr.pl'
vcf_basic = '/users/selasady/my_titan_itmi/impala_scripts/annotation/parse_vcf.pl'

# load packages
from impala.dbapi import connect
from impala.util import as_pandas
import time
import csv
import subprocess
import numpy as np
import sys
import datetime
now = datetime.datetime.now()

# TODO: update from test_vars to global_vars once table is complete
# create query to download variants from input table that are not in global_vars
comparison_query = '''
with vars as (
  select chrom, pos, id as rs_id, ref, alt
  from p7_platform.wgs_illumina_variant
  where chrom = '1'
  and id is not null)

select chrom, pos, rs_id, ref, alt
from vars t
where not exists (
  select chrom, pos, rs_id, ref, alt
  from p7_product.all_variants g
  where g.chrom = '1'
  and t.pos = g.pos
  and t.ref = g.ref
  and t.alt = g.alt)
limit 15
'''
# select chrom, pos, rs_id, ref, alt
# from {}.{} t
# where not exists (
#   select chrom, pos, rs_id, ref, alt
#   from p7_product.test_vars g
#   where t.chrom = g.chrom
#   and t.pos = g.pos
#   and t.ref = g.ref
#   and t.alt = g.alt
#   )'''.format(input_db, input_table)

#connect to impala
conn=connect(host='glados19', port=21050)
#create a cursor object to interact with db
cur = conn.cursor()
# execute sql query
print "Searching for variants that are not in the global variants table... \n"
cur.execute(comparison_query)
# save results as pandas dataframe
new_vars = as_pandas(cur)


#################################################
## create vcf files by row for each chromosome ##
#################################################
# create vcf header
def create_header(outfile_name):
    print "Creating VCF header. \n"
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

# function to create vcf file
def create_vcf(out_name):
    print "Creating VCF files. \n"
    # create named vcf file
    vcf_out = out_name + '.vcf'
    # create header for file
    create_header(vcf_out)
    # write variants to file row by row to save memory
    try:
        new_vars.to_csv(vcf_out, sep='\t', encoding='utf-8', mode='a', header=False, index=False)
    except Exception as csv_error:
        print csv_error

# function to verify vcf format using GATK's barebones.pl script
def check_vcf(out_name):
    print "Verifying VCF format. \n"
    vcf_in =  out_name + '.vcf'
    vcf_out = out_name + '_verified.vcf'
    # create the file and run snpeff
    with open(vcf_out, "w") as out_file:
        try:
            subprocess.call(['perl', vcf_basic, vcf_in], stdout=out_file)
        except subprocess.CalledProcessError as e:
             print e.output

# function to run verified vcf files through snpeff
def run_snpeff(out_name):
    print "Annotating variants with coding consequences using snpeff. \n"
    vcf_in = out_name + '_verified.vcf'
    vcf_out = out_name + '_snpeff.vcf'
    with open(vcf_out, "w") as f:
        try:
            subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "-t", "-v", "-noStats", "GRCh37.75", vcf_in], stdout=f)
        except subprocess.CalledProcessError as e:
             print e.output

# output snpeff effects as tsv file with one effect per line
def parse_snpeff(out_name):
    print "Parsing snpeff output. \n"
    vcf_in = out_name + '_snpeff.vcf'
    tsv_out = out_name + '.tsv'
    # create command to parse snpeff
    snpout_cmd = 'cat {} | {} | {} -jar {} extractFields \
    - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
    "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
    "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(vcf_in, snpeff_oneperline_perl, \
    java_path, snpsift_jar, tsv_out)
    # call subprocess and communicate to pipe output between commands
    ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    print ps.communicate()[0]

# remove header so strange characters don't break impala
def remove_header(out_name):
    print "Removing header for upload to impala. \n"
    tsv_in = out_name + '.tsv'
    tsv_out = out_name + '_final.tsv'
    tsv_cmd = "sed '1d' {} > {}".format(tsv_in,tsv_out)
    tsv_proc = subprocess.Popen(tsv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print tsv_proc.communicate()[0]

# upload results to hdfs
def upload_hdfs(out_name):
    # define output path on hdfs
    out_path = "{}snpeff_{}".format(hdfs_path, str(now.strftime("%Y%m%d")))
    mkdir_cmd = "hdfs dfs -mkdir {}".format(out_path)
    mkdir_proc = subprocess.Popen(mkdir_cmd, shell=True, stderr=subprocess.STDOUT)
    print mkdir_proc.communicate()[0]
    # put snpeff file in the newly created directory
    print "Uploading coding consequences to HDFS... \n"
    snpeff_file = out_name + '_final.tsv'
    hdfs_cmd = 'hdfs dfs -put {} {}'.format(snpeff_file, out_path)
    hdfs_proc = subprocess.Popen(hdfs_cmd, shell=True, stderr=subprocess.STDOUT)
    print hdfs_proc.communicate()[0]
    # set read/write permissions on directory
    chown_dir_cmd = "hdfs dfs -chown -R impala:supergroup {}".format(hdfs_path)
    chown_proc = subprocess.Popen(chown_dir_cmd, shell=True, stderr=subprocess.STDOUT)
    print chown_proc.communicate()[0]

def create_table(out_name):
    print "Creating table to store results. \n"
    table_name = out_name + '_' + str(now.strftime("%Y%m%d"))
    # create empty table to insert annotated new variants
    create_coding_table = '''
        create table {}.{}
         (chrom string,
          pos int,
          rs_id string,
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
    '''.format(input_db, table_name)
    try:
        cur.execute(create_coding_table)
    except subprocess.CalledProcessError as e:
             print e.output

# load hdfs files into table
def results_to_table(out_name):
    print "Loading results to impala table. \n"
    out_path = "{}snpeff_{}".format(hdfs_path, str(now.strftime("%Y%m%d")))
    table_name = out_name + '_' + str(now.strftime("%Y%m%d"))
    load_query = '''
        load data inpath '{}' into table {}.{}
    '''.format(out_path, input_db, table_name)
    try:
        cur.execute(load_query)
    except subprocess.CalledProcessError as e:
             print e.output

# compute stats on newly created table
def stats_coding(out_name):
    table_name = out_name + '_' + str(now.strftime("%Y%m%d"))
    print "Running compute stats on {}".format(table_name)
    coding_compstats = "compute stats {}.{}".format(input_db, table_name)
    try:
        cur.execute(coding_compstats)
    except subprocess.CalledProcessError as e:
             print e.output

# annotate variants with ensembl

# if new variants are found, annotate with snpeff and upload to impala as a table
if len(new_vars) > 0:
    print str(len(result_name)) + " new variant(s) were found. \n"
    create_vcf(result_name)
    check_vcf(result_name)
    run_snpeff(result_name)
    parse_snpeff(result_name)
    # remove_header(result_name)
    # upload_hdfs(result_name)
    # create_table(result_name)
    # results_to_table(result_name)
    # stats_coding(result_name)

    sys.exit("New variants added to global variants table.")
    cur.close()

else:
    sys.exit("No new variants found.")
    cur.close()

# annoate with clinvar





# run through snpeff
# add to table