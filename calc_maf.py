__author__ = 'selasady'

#######################
# connect to database #
#######################
print "Connecting to impala..." + "\n"
from impala.dbapi import connect

conn = connect(host='glados19', port=21050)
cur = conn.cursor()

#test query
query = "SELECT cgi.sample_id, concat(cgi.chr, ':', CAST(cgi.start AS STRING), ':', CAST(cgi.stop AS STRING)) as pos, cgi.ref, cgi.allele1seq, cgi.allele2seq, CASE WHEN cgi.allele1seq <> cgi.ref THEN '1'  ELSE '0' END as allele1_MA_count, CASE WHEN cgi.allele2seq <> cgi.ref THEN '1'  ELSE '0' END as allele2_MA_count FROM p7_ptb.comgen_variant cgi LIMIT 500"

#########################
# Run Query on Impala ##
#########################
print "Running query on impala..." + "\n"
from impala.util import as_pandas
cur.execute(query)
query_results = as_pandas(cur)
if len(query_results) > 0:
    print "Here's a preview of the results: \n"
    print query_results.head(5)
else:
    print "No results found"


# 0 = ref
# 1 = alt
#count total number of alt at that position and divide by
# total alleles = 2 alleles * n_samples

#############
# Calc MAF ##
#############
tot_alleles =  len(query_results) * 2
