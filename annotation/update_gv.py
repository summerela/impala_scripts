# specify table to check for new variants
input_db = "p7_product"
input_table =  "test"

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

# TODO: update from global_vars to global_variants once table is complete
# create query to download variants from input table that are not in global_vars
comparison_query = '''
select chrom, pos, ref, alt
from {}.{} t
where not exists (
  select chrom, pos, ref, alt
  from p7_product.global_vars g
  where t.chrom = g.chrom
  and t.pos = g.pos
  and t.ref = g.ref
  and t.alt = g.alt
  )'''.format(input_db, input_table)

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
    # create named vcf file
    vcf_out = out_name + '.vcf'
    # create header for file
    create_header(vcf_out)
    # write variants to file row by row to save memory
    new_vars.to_csv(vcf_out, sep='\t', encoding='utf-8', mode='a')

# function to verify vcf format using GATK's barebones.pl script
def check_vcf(out_name):
    vcf_in =  out_name + '.vcf'
    vcf_out = out_name + '_verified.vcf'
    # create the file and run snpeff
    with open(vcf_out, "w") as out_file:
        try:
            subprocess.call(['perl', vcf_basic, vcf_in], stdout=out_file)
        except subprocess.CalledProcessError as e:
             print e.output

# if new variants are found, annotate with snpeff and upload to impala as a table
if len(new_vars) > 0:
    print str(len(new_vars)) + " new variant(s) were found. \n"
    print "Creating VCF files. \n"
    try:
        create_vcf("new_vars")
    except Exception as create_error:
        print create_error
    print "Verifying VCF format. \n"
    try:
        check_vcf("new_vars")
    except Exception as check_error:
        print check_error
    sys.exit("New variants added to global variants table.")
else:
    sys.exit("No new variants found.")

# annoate with clinvar





# run through snpeff
# add to table