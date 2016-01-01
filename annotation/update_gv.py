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
chroms = map( str, range(1,23) ) + ['X','Y','MT', 'M']

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
    create_header(vcf_out)
    # connect to vars_to_snpeff table
    # write variants to file row by row to save memory
    with open(vcf_out, 'a') as csvfile:
        for row in new_vars:
            writer = csv.writer(csvfile, delimiter="\t", lineterminator = '\n')
            writer.writerow(row)

# function to create vcf file
def create_vcf(db_name, table_name, chrom_name, out_name):
    # create named file for each chromosome
    vcf_out = 'chr' + chrom_name + '_' + out_name + '.vcf'
    # create header for each chromosome file
    # TODO: make this one function instead of calling header function
    create_header(vcf_out)
    # connect to vars_to_snpeff table
    # write variants to file row by row to save memory
    with open(vcf_out, 'a') as csvfile:
        for row in new_vars:
            writer = csv.writer(csvfile, delimiter="\t", lineterminator = '\n')
            writer.writerow(row)

# function to verify vcf format using GATK's barebones.pl script
def check_vcf(chrom_name, out_name):
    vcf_in =  'chr' + chrom_name + '_' + out_name + '.vcf'
    vcf_out = 'chr' + chrom_name + '_verified.vcf'
    # create the file and run snpeff
    with open(vcf_out, "w") as out_file:
        try:
            subprocess.call(['perl', vcf_basic, vcf_in], stdout=out_file)
        except subprocess.CalledProcessError as e:
             print e.output

# if new variants are found, annotate with snpeff and upload to impala as a table
if len(new_vars) > 0:
    print str(len(new_vars)) + " new variant(s) were found. \n"
    for chrom in chroms:
        print "Creating VCF files for chromosome {}. \n".format(chrom)
        create_vcf(input_db, input_table, chrom, "new_vars")
        print "Verifying VCF format for chromosome {}. \n".format(chrom)
        check_vcf(chrom, "new_vars")
else:
    print "No new variants found."

# annoate with clinvar





# run through snpeff
# add to table