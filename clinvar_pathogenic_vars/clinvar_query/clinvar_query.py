__author__ = 'Summer Elasady for ISB'

#################
# load modules  #
#################
import argparse
import sys
import pandas as pd
from impala.dbapi import connect
import re
from impala.util import as_pandas

#ignore irrelevant pandas warnings
pd.options.mode.chained_assignment = None

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
    #####################
    ## Generate Report ##
    #####################
    print "Generating an html report..."
    # make sql query pretty
    rep = {'FROM':'<br> FROM', 'JOIN':'<br>JOIN', 'ON':'<br>ON', 'AND':'<br>AND', 'WHERE':'<br>WHERE'}
    rep = dict((re.escape(k), v) for k, v in rep.iteritems())
    pattern = re.compile("|".join(rep.keys()))
    pretty_query= pattern.sub(lambda m: rep[re.escape(m.group(0))], query)
    #make results pretty
    summary_table_1 = query_df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">')
    #Create html report
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <h1>Clinvar Annotation of Variants</h1>
            <h2>Specifications: </h2>
                <p> Clinvar was used to annotate variants falling in regions marked with a clinical significance of 4 or 5,
                but not 2 or 3. The following paramaters were specified: <br>
                <ul>
                <li>Chromosome(s) of interest =  ''' + args.chr + '''</li>
                <li>Gene(s) of interest =  ''' + args.genes + '''</li>
                <li>Zygosity =  ''' + args.zygosity + '''</li>
                <li>Trio Member(s) =  ''' + args.member + '''</li>
                <li>Sample ID'(s) =  ''' + args.sample_id + '''</li>
                <li>Platform =  ''' + args.platform + '''</li>
                </ul>
            <!-- *** Query *** --->
                <h2>Query: </h2>
                <p>The following query was run on impala: <br><br>''' + pretty_query + '''</p>
            <!-- *** Results *** --->
                <h2>Results: </h2>
                <p>The following results were found:<br>
                <ul>
                <li>Total Variants Found:''' + len(query_df) + '''</li>
                <li>Total Variants Found:''' + pd.unique(df.column_name.ravel()) + '''</li>
                </p>
        </body>
        </html>'''
else:
    print "No results found"




f = open('test.html','w')
f.write(html_string)
f.close()

#close all the things
cur.close()
sys.exit()



### TO DO ###
# illumina variants file needs gt,zygosity fields added
# also need to add member column to both cgi and illumina in impala
# parse zygosity depending on platform
# add filtering options
# link to gene info
# link to dbSNP
# scrollable, filterable table