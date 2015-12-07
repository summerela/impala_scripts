__author__ = 'Summer Elasady for ISB'

#################
# load modules  #
#################
import argparse
import sys
import pandas as pd
from impala.dbapi import connect
from impala.util import as_pandas

#ignore irrelevant pandas warnings
pd.options.mode.chained_assignment = None

###################
# parse user args #
###################
try:
    parser = argparse.ArgumentParser(description='Find pathogenic variants using clinvar.', usage='%(prog)s [options]')
    parser.add_argument("--chr", help="chromosome as 7,8,9 etc. or 'all'", type=str, nargs='?', default='all')
    parser.add_argument("--zygosity", help="zygosity as hom,het-ref,het-alt or 'all'", type=str, nargs='?', default='all')
    parser.add_argument("--genes", help="genes as HOXB1,BRCA1,BRCA2 etc. or 'all'", type=str, nargs='?', default='all')
    parser.add_argument("--member", help="trio members as M,F,NB or 'all'", type=str, nargs='?', default='all')
    parser.add_argument("--sample_id", help="sample id's as 101-454-M,101-345-F etc. or 'all'", type=str, nargs='?',
                        default='all')
    parser.add_argument("--platform", help="cgi or illumina", type=str, nargs='?', default='illumina')
    args = parser.parse_args()
except:
    e = sys.exc_info()[0]
    print e

####################
# format user args #
####################
print "Creating a query with the arguments you provided..."

# create empty list to store query conditions
conditionals = []

# function to process user args
def process_args(arg, val):
    # if arg is not 'all' and there are more than one argument
    if val != 'all' and (',' in val):
        conditionals.append(" vcf.{0} IN ('".format(arg) + "','".join(val.split(',')) + "')")
    # if arg is not 'all' but there is only one arg
    elif val!= 'all' and  (',' not in val):
        conditionals.append(" vcf.{0} = '{1}'".format(arg,val))
    #if arg is 'all' leave blank
    elif val == 'all':
        pass
    #if arg not any of above, error and halt
    else:
        sys.exit("Check that your command line args are properly formatted and try again.")

# run process_args on user args
for key, value in vars(args).items():
    #dont add the platform argument to the query
    if key != 'platform':
        process_args(key, value)

# if there is more than one user argument
if len(conditionals) > 1:
    #first argument gets 'WHERE' and the rest get 'AND'
    query_conditions = ['AND' + s for s in conditionals]
    query_args = " ".join(map(str, query_conditions))
# if there is only one user arg, it gets surrounded by quotes
else:
    query_args = ['AND' + conditionals]

###############
# Build Query #
###############
print "Building a query to look for " + args.zygosity + " " + args.platform + " variants in gene(s):" + str(args.genes) + ", for member(s):" + args.member+\
    " in sample(s):" + args.sample_id + " on chr:" + args.chr+ "."

#build query
if args.platform == 'cgi':
    query = """
SELECT vcf.chr, vcf.start, vcf.stop, vcf.ref, vcf.allele1seq as alt1, vcf.allele2seq as alt2, vcf.zygosity, vcf.vartype,
vcf.allele1varquality as alt1_qual, vcf.allele2varquality as alt2_qual, vcf.totalreadcount, vcf.sample_id,
clin.rsid, clin.rspos, clin.geneinfo, clin.dbsnpbuildid as dbsnp_build, clin.hgvs, clin.clin_sig
    FROM public_hg19.clinvar_grch37 AS clin
    JOIN p7_platform.wgs_comgen_variant AS vcf
    ON clin.chr = vcf.chr
    AND clin.ref = vcf.ref
    AND clin.pos BETWEEN vcf.start and vcf.stop
    AND (clin.alt = vcf.allele1seq or clin.alt = vcf.allele2seq)
WHERE ((FIND_IN_SET('4', clin.clin_sig) > 0) OR (FIND_IN_SET('5', clin.clin_sig)>0)
    AND
       (FIND_IN_SET('3', clin.clin_sig)=0 AND FIND_IN_SET('2', clin.clin_sig)=0))
{0:s}
             """.format(query_args)
elif args.platform == 'illumina':
    query = """
SELECT vcf.chr, vcf.pos, vcf.ref, vcf.alt, vcf.qual, vcf.filter, vcf.gt, vcf.zygosity, vcf.sample_id,
clin.rsid, clin.rspos, clin.geneinfo, clin.dbsnpbuildid as dbsnp_build, clin.hgvs, clin.clin_sig
 FROM public_hg19.clinvar_grch37 AS clin
 JOIN p7_platform.wgs_illumina_variant AS vcf
    ON clin.chr = vcf.chr
    AND clin.ref = vcf.ref
    AND clin.pos = vcf.pos
    AND clin.alt = vcf.alt
WHERE ((FIND_IN_SET('4', clin.clin_sig) > 0) OR (FIND_IN_SET('5', clin.clin_sig)>0)
    AND
       (FIND_IN_SET('3', clin.clin_sig)=0 AND FIND_IN_SET('2', clin.clin_sig)=0))
{0:s}
            """.format(query_args)
else:
    print "Did you select illumina or cgi as your platform? Please check and try again."
print query
######################
# connect to database #
######################
print "Connecting to impala..." + "\n"

#create database connection
conn = connect(host='glados19', port=21050)
cur = conn.cursor()

#########################
# Run Query on Impala ##
#########################
print "Running the following query on impala..." + "\n"
print query

#execute query on impala
cur.execute(query)
query_df = as_pandas(cur)

#check that results are returned and print preview
if len(query_df) > 0:
    print "Results saved to current working directory as clinvar_results.csv..."
    #save to csv
    query_df.to_csv('clinvar_results.csv', header=True, encoding='utf-8', index=False)
else:
    print "No results found"




#close all the things
cur.close()
sys.exit()



### TO DO ###
# illumina variants file needs gt,zygosity fields added
# also need to add member column to both cgi and illumina in impala
# parse zygosity depending on platform