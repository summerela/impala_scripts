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

# declare script variable
chroms = map( str, range(1,23) ) + ['X','Y','MT', 'M']

#connect to impala
conn=connect(host='glados16', port=21050)
#create a cursor object to interact with db
cur = conn.cursor()

##########################################
## search for variants not in gv table  ##
##########################################
def get_vars(input_db, input_table, chrom_name):
# create query to download variants from input table that are not in global_vars
    comparison_query = '''
    select chrom, pos, rs_id, ref, alt
    from {}.{} t
    where chrom = '{}'
    and not exists (
      select chrom, pos, rs_id, ref, alt
      from p7_product.test_vars g
      where t.chrom = g.chrom
      and t.pos = g.pos
      and t.ref = g.ref
      and t.alt = g.alt
      )
      order by pos
      '''.format(input_db, input_table, chrom_name)
    # execute sql query
    cur.execute(comparison_query)
    # store results as pandas table
    results = as_pandas(cur)
    return results

##########################################
## create vcf files for each chromosome ##
##########################################
# create vcf header
def create_header(outfile_name):
   # create vcf header
    lines=[]
    lines.append('##fileformat=VCFv4.0')
    lines.append('##fileDate='+ time.strftime("%y%m%d"))
    lines.append('##reference=grch37 v.75')
    lines.append('#CHROM\t' + 'POS\t' + 'ID\t' + 'REF\t' + 'ALT\t' + 'QUAL\t'+ 'FILTER\t' + 'INFO\t' + 'FORMAT\t' + 'SAMPLE\t' + '\n')
    header = '\n'.join(lines)
    out = open(outfile_name, 'wb')
    out.write(header)
    out.close()

# function to create vcf file
def create_vcf(out_name, chrom, var_df):
    # create named vcf file
    vcf_out = "chr" + str(chrom) + '_' + out_name + '.vcf'
    print "\n Creating VCF files for chromosome {}... \n".format(chrom)
    # create header for file
    create_header(vcf_out)
    # write variants to file
    try:
        var_df.to_csv(vcf_out, mode='a', sep='\t', header=False, index=False)
    except error as e:
        print e

# function to verify vcf format using GATK's barebones.pl script
def check_vcf(out_name):
    print "Verifying VCF format for chromosome {}. \n".format(chrom)
    vcf_in =  "chr" + str(chrom) + '_' + out_name + '.vcf'
    vcf_out = "chr" + str(chrom) + '_' + out_name + '_verified.vcf'
    # create the file and run snpeff
    with open(vcf_out, "w") as out_file:
        try:
            subprocess.call(['perl', vcf_basic, vcf_in], stdout=out_file)
        except subprocess.CalledProcessError as e:
             print e.output

# function to run verified vcf files through snpeff
def run_snpeff(out_name):
    print "Adding coding consequences for chromosome {}. \n".format(chrom)
    vcf_in = "chr" + str(chrom) + '_' + out_name + '_verified.vcf'
    vcf_out = "chr" + str(chrom) + '_' + out_name + '_snpeff.vcf'
    with open(vcf_out, "w") as f:
        try:
            subprocess.call([java_path, "-Xmx16g", "-jar", snpeff_jar, "-t", "-v", "-noStats", "GRCh37.75", vcf_in], stdout=f)
        except subprocess.CalledProcessError as e:
             print e.output

# output snpeff effects as tsv file with one effect per line
def parse_snpeff(out_name):
    print "Parsing snpeff output for chromosome {}. \n".format(chrom)
    vcf_in = "chr" + str(chrom) + '_' + out_name + '_snpeff.vcf'
    tsv_out = "chr" + str(chrom) + '_' + out_name + '.tsv'
    # create command to parse snpeff
    snpout_cmd = 'cat {} | {} | {} -jar {} extractFields \
    - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
    "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
    "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {}'.format(vcf_in, snpeff_oneperline_perl, \
    java_path, snpsift_jar, tsv_out)
    # call subprocess and communicate to pipe output between commands
    try:
        ps = subprocess.Popen(snpout_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        print ps.communicate()[0]
    except:
        e = sys.exc_info()[0]
        print e

# remove header so strange characters don't break impala
def remove_header(out_name):
    print "Removing header on chromosome {} for upload to impala. \n".format(chrom)
    tsv_in = "chr" + str(chrom) + '_' + out_name + '.tsv'
    tsv_out = "chr" + str(chrom) + '_' + out_name + '_final.tsv'
    tsv_cmd = "sed '1d' {} > {}".format(tsv_in,tsv_out)
    try:
        tsv_proc = subprocess.Popen(tsv_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print tsv_proc.communicate()[0]
    except Exception as e:
        print e

# process new variants by chromosome
for chrom in chroms:
    new_vars = get_vars(input_db, input_table, chrom)
    if len(new_vars) > 0:
        create_vcf(result_name, chrom, new_vars)
        check_vcf(result_name)
        run_snpeff(result_name)
        parse_snpeff(result_name)
        remove_header(result_name)



cur.close()






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
    except Exception as e:
             print e

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
    except Exception as e:
             print e

# compute stats on newly created table
def stats_coding(out_name):
    table_name = out_name + '_' + str(now.strftime("%Y%m%d"))
    print "Running compute stats on {}".format(table_name)
    coding_compstats = "compute stats {}.{}".format(input_db, table_name)
    try:
        cur.execute(coding_compstats)
    except Exception as e:
             print e

# annotate variants with ensembl




# annoate with kaviar





# run through snpeff
# add to table