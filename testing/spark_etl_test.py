#!/usr/bin/env pyspark

'''
Assumptions:
Input file is in impala table at least the following columns:
- chrom
- pos
- ref
- alt
- subject_id
- sample_id

To Run:

- Edit main routine

Type the following from command line

pyspark spark_etl_test.py  2>/dev/null


'''

import os
import unittest

from pyspark import SparkContext, SparkConf, SQLContext


############################
### Unit Test with Spark ###
############################

class TestETL(unittest.TestCase):

    def __init__(self, in_table, local_dir='./', hdfs_dir='/vlad/wgs_ilmn.db', appname='spark_app'):
        self.in_table = in_table
        self.local_dir = local_dir
        self.hdfs_dir = hdfs_dir
        self.appname = appname
        self.conf = (SparkConf()
               .setAppName(self.appname))
        self.sc = SparkContext(conf=self.conf)
        self.somatic_chr = sorted(map(str, range(1,23) + ["M", "X", "Y"]))
        self.sqlContext = SQLContext(self.sc)
        # encode any parquet binary formatting as strings
        self.sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")
        # # cache table metadata for query speed
        self.sqlContext.sql("SET spark.sql.parquet.cacheMetadata=true")
        self.data = self.sqlContext.read.parquet("{}".format(self.in_table))
        # register parquet file as temp table to query
        self.parquet_tbl = self.data.registerTempTable("parquet_tbl")


    def check_empty(self, col_name):
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE {} IS NULL".format(col_name)).count()
        self.assertTrue(result == 0, "Null ?values found in {} column.".format(col_name))

    def check_chrom_set(self):
        result = self.sqlContext.sql("SELECT distinct chrom FROM parquet_tbl")
        uploaded_chroms = sorted(map(str, result.rdd.map(lambda r: r.chrom).collect()))
        self.assertListEqual(self.somatic_chr, uploaded_chroms,
                             "The following chromosomes are not loaded: \n {}".format(
                                 set(self.somatic_chr).difference(uploaded_chroms)))

    def check_chrom_count(self):
        result = self.sqlContext.sql("SELECT chrom, COUNT(*) FROM parquet_tbl GROUP BY chrom HAVING COUNT(*) < 10000")
        self.assertTrue(result.count() == 0, "Not all chromosomes contain at least 10000 rows.")

    def check_chrom_prefix(self):
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE chrom LIKE '*chr*'")
        self.assertTrue(result.count() == 0, "The 'chr' prefix was not removed from the chrom column.")

    def check_chrom_mt(self):
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE chrom = 'MT' or chrom = 'mt' ")
        self.assertTrue(result.count() == 0, "The 'MT' prefix was not changed to 'M' in the chrom column.")

    def ref_check(self):
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE regexp_extract( ref, '([^ACGTN.])', 0 ) IS NULL")
        self.assertTrue(result.count() == 0, "The ref field contains values other than A,T,G,C or N.")

    def test_alleles(self):
        result = self.sqlContext.sql("SELECT chrom, pos, ref, alt, subject_id, sample_id FROM parquet_tbl WHERE \
            length(ref) = length(alt) AND alt <> '.' GROUP BY chrom, pos, ref, alt, subject_id, sample_id \
        HAVING COUNT(*) > 2")
        self.assertTrue(result.count() == 0, "The alt field contains values other than A,T,G,C or N.")

    def test_ref(self):
        result = self.sqlContext.sql("SELECT chrom, pos, ref FROM parquet_tbl WHERE length(ref) = length(alt) \
            AND ref <> '.' GROUP BY chrom, pos, ref HAVING COUNT(*) > 2")
        # self.assertTrue(result.count() == 0, "Multiple values for the reference field were found at: {}".format(result))
        print result.collect()


    def tear_down(self):
        # clear memory cache
        self.sqlContext.clearCache()
        # close connection
        self.sc.stop()

    def run_tests(self):

        try:
            os.path.isfile(self.in_table)
        except Exception as e:
            print e




        # check chrom, pos, sample_id column vals not empty
        # self.check_empty("chrom")
        # self.check_empty("ref")
        # self.check_empty("subject_id")
        # self.check_empty("sample_id")

        # check that all somatic chroms loaded
        # self.check_chrom_set()

        # # check that each chrom contains at least 1000 rows
        # self.check_chrom_count()

        # # check that 'chr' prefix has been removed from chromosome
        # self.check_chrom_prefix()
        #
        # # check that mitochondrial chromosome has been renamed from MT to M
        # self.check_chrom_mt()
        #
        # # alt cols contains expected values A T G C N only
        # self.ref_check()
        #
        # # # check for no more than two alleles per person at chrom, pos, ref
        # self.test_alleles()
        # #
        # # # # ref allele should be the same at each position
        self.test_ref()

        self.tear_down()


###############
if __name__ == '__main__':

    # instantiate class
    spark = TestETL(in_table='hdfs://nameservice1/vlad/wgs_ilmn.db/vcf_variant/', local_dir=os.getcwd(), hdfs_dir='/user/selasady/', appname='vcf_testing')


    # run tests
    spark.run_tests()




######################## unused code snippets #######################

# def tsv_to_df(self, input_file):
#     # import file as dataframe, all cols will be imported as strings
#     df = self.sqlContext.read.format("com.databricks.spark.csv").option("header", "true").option("delimiter", "\t").option("inferSchema", "true").load(input_file)
#     # # cache df object to avoid rebuilding each time
#     df.cache()
#     # register as temp table for querying
#     df.registerTempTable("spark_df")
#     return df

# def get_ref_snp(self, input_df):
#     random_snp = input_df[['chrom', 'pos', 'ref', 'alt']][input_df['ref'] != '.'].sample(1)
#     return random_snp
#
# def query_snp(self, input_snp):
#     # testing if variant has 1-based coords
#     exclude = ['ins', 'del']
#     query_string = "chr{}:{}".format(input_snp.iloc[0]['chrom'], input_snp.iloc[0]['pos'])
#     test_url = 'http://myvariant.info/v1/query?q={}'.format(query_string)
#     response = urllib.urlopen(test_url)
#     data = json.loads(response.read())
#     ref_var = []
#     while len(ref_var) == 0:
#         for x in range(0,4):
#             try:
#                  ref_hit = data['hits'][0]['_id'] if len(data['hits'][0]['_id']) > 0 else 'null'
#                  if any(x in exclude for x in ref_hit):
#                      pass
#             except:
#                  time.sleep(10)
#         else:
#             ref_var.append(ref_hit)
#             ref_out = ''.join(ref_var)
#             return ref_out
#
# def check_one_based(self, input_df):
#     test_snp = self.get_ref_snp(input_df)
#     ref_snp = self.query_snp(test_snp)
#     file_snp =  "chr{}:g.{}{}>{}".format(test_snp.iloc[0]['chrom'], test_snp.iloc[0]['pos'],  test_snp.iloc[0]['ref'],  test_snp.iloc[0]['alt'])
#     self.assertEqual(ref_snp, file_snp, "File variant {} does not match reference {}".format(file_snp, ref_snp))
