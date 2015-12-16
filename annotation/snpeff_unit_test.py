# impala connection strings
impala_host = 'glados18'
impala_port = '21050'

# Database to test
input_db = "p7_product"
input_table = "all_coding"

####################
## import modules ##
####################
import pandas as pd
from impala.dbapi import connect
from impala.util import as_pandas
import numpy as np
import urllib, json

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

## create connection to impala
conn=connect(host=impala_host, port=impala_port, timeout=10000, user='selasady')
cur = conn.cursor()

#################
# Unit Testing ##
#################

# define list of possible expected chromsomes, mitochondria are not always present
chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', \
              '21', '22', 'X', 'Y', 'M', 'MT']

def get_stats(in_db, in_table):
    # create object to get table column stats for use with downstream functions
    get_col_stats = "show column stats {}.{}".format(in_db, in_table)
    cur.execute(get_col_stats)
    col_stats_df = as_pandas(cur)
    col_stats_df = col_stats_df.replace('-1', 'NA')
    col_stats_df.columns = ['column', 'data_type', 'distinct_values', 'num_nulls', 'max_size', 'avg_size']
    return col_stats_df

# create variables to test
stats_df = get_stats(input_db, input_table)
num_chroms = stats_df['distinct_values'][stats_df['column'] == 'chrom'].iloc[0]
chrom_dtype = stats_df['data_type'][stats_df['column'] == 'chrom'].iloc[0]
pos_dtype = stats_df['data_type'][stats_df['column'] == 'pos'].iloc[0]

# check if any columns contain only a single value or all-nulls
non_distincts =  stats_df.loc[stats_df['distinct_values'] < 2]
if len(non_distincts) > 0:
    print "FAIL: The following columns do not have more than one distcint value: \n"
    print non_distincts
else:
    print "PASS: All columns have more than one distinct value. \n"

# check that chromosome column contains 24 chromosomes
if num_chroms != 24:
    print "FAIL: " + str(num_chroms) + " chromosomes found. \n"
    print sorted(stats_df['chrom'])
else:
    print "PASS: Table contains 24 chromosomes. \n"

# check that chromosome column is a string
if chrom_dtype != 'STRING':
    print "FAIL: Chromosome column is listed as {} not a string. \n".format(chrom_dtype)
else:
    print "PASS: Chromosome data type is a string, as expected. \n"

# check that pos column is an integer
if pos_dtype != 'INT':
    print "FAIL: Pos column is listed as {} not an integer. \n".format(pos_dtype)
else:
    print "PASS: Position column was an integer, as expected. \n"

# test that all chromosomes are represented
def test_chroms_present(in_db, in_table):
    chrom_query = "select distinct chrom from {}.{}".format(in_db, in_table)
    cur.execute(chrom_query)
    chrom_list = sorted(cur.fetchall())
    chrom_diffs = np.setdiff1d(np.array(chroms),np.array(chrom_list))
    if len(chrom_diffs) > 0:
        print "FAIL: The following chromosomes were not found: {}. \n".format(chrom_diffs)
    else:
        print "PASS: All chromosomes were successfully uploaded. \n"

test_chroms_present(input_db, input_table )

# test that coordinates are 1-based and match online source
def test_against_online(in_db, in_table):
    # get random row from table
    grab_row = "select * from {}.{} order by rand () limit 1".format(in_db, in_table)
    cur.execute(grab_row)
    tester_row = as_pandas(cur)
    # use rsid to get information from web api
    test_url = 'http://myvariant.info/v1/query?q={}'.format(tester_row['id'])
    response = urllib.urlopen(test_url)
    data = json.loads(response.read())
    # test that chromosome match
    ref_chrom = data["hits"][0]["dbsnp"]["chrom"]
    if str(tester_row['chrom'].iloc[0]) == str(ref_chrom):
        print "PASS: chromosome matches online reference for {}. \n".format(tester_row['id'].iloc[0])
    else:
        print "FAIL: chromosome does not match online reference for {}. \n".format(tester_row['id'].iloc[0])
    #test that position matches online reference
    ref_position =  data["hits"][0]["vcf"]["position"]
    if str(tester_row['pos'].iloc[0]) == str(ref_position):
        print "PASS: position matches online reference for {}. \n".format(tester_row['id'].iloc[0])
    else:
        print "FAIL: position does not match 1-based online reference for {}. \n".format(tester_row['id'].iloc[0])
    #test that ref allele matches online reference allele
    ref_ref =  data["hits"][0]["vcf"]["ref"]
    if str(tester_row['ref'].iloc[0]) == str(ref_ref):
        print "PASS: reference allele matches online reference for  {}. \n".format(tester_row['id'].iloc[0])
    else:
        print "FAIL: reference allele does not match online reference for  {}. \n".format(tester_row['id'].iloc[0])
    # test that alt allele matches alternate allele of online reference
    ref_alt = position =  data["hits"][0]["vcf"]["alt"]
    if str(tester_row['alt'].iloc[0]) == str(ref_alt):
        print "PASS: alt allele matches online reference for  {}. \n".format(tester_row['id'].iloc[0])
    else:
        print "FAIL: alt allele does not match online reference for  {}. \n".format(tester_row['id'].iloc[0])
    #make sure that gene name matches snpeff gene name for this rsid
    ref_gene_name = data["hits"][0]["snpeff"]["ann"][0]["gene_id"]
    if tester_row['gene'].iloc[0] == ref_gene_name:
         print "PASS: snpeff gene name {} matches online reference for {}. \n".format(ref_gene_name, tester_row['id'].iloc[0])
    else:
       print "FAIL snpeff gene_name {} does not match online reference for {}. \n".format(ref_gene_name, tester_row['id'].iloc[0])

test_against_online(input_db, input_table)

cur.close()