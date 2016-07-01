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

    def __init__(self, in_table, is_vcf='no', platform='illumina', local_dir='./', hdfs_dir='/vlad/wgs_ilmn.db', appname='spark_app'):
        self.in_table = in_table
        self.is_vcf = is_vcf
        self.platform = platform
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
        # check for null values
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE {} IS NULL".format(col_name)).count()
        self.assertTrue(result == 0, "Null ?values found in {} column.".format(col_name))

    def check_chrom_set(self):
        # check that all chroms are loaded
        result = self.sqlContext.sql("SELECT distinct chrom FROM parquet_tbl")
        uploaded_chroms = sorted(map(str, result.rdd.map(lambda r: r.chrom).collect()))
        self.assertListEqual(self.somatic_chr, uploaded_chroms,
                             "The following chromosomes are not loaded: \n {}".format(
                                 set(self.somatic_chr).difference(uploaded_chroms)))

    def check_chrom_count(self):
        # check that data was uploaded for each chrom
        result = self.sqlContext.sql("SELECT chrom, COUNT(*) FROM parquet_tbl GROUP BY chrom HAVING COUNT(*) < 10000")
        self.assertTrue(result.count() == 0, "Not all chromosomes contain at least 10000 rows.")

    def check_chrom_prefix(self):
        # check that 'chr' prefix was removed in chrom column
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE chrom LIKE '*chr*'")
        self.assertTrue(result.count() == 0, "The 'chr' prefix was not removed from the chrom column.")

    def check_chrom_mt(self):
        # check that mitochondria was renamed from MT to M
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE chrom = 'MT' or chrom = 'mt' ")
        self.assertTrue(result.count() == 0, "The 'MT' prefix was not changed to 'M' in the chrom column.")

    def allele_check(self):
        # check that the allele column only contains A,T,G,C, N or .
        result = self.sqlContext.sql("SELECT allele, regexp_replace(allele, 'A|T|G|C|N|\\.', '') as test_col \
        from parquet_tbl WHERE regexp_replace(allele, 'A|T|G|C|N|\\.', '')  <> ''")
        self.assertTrue(result.count() == 0, "The ref field contains characters other than A,T,G,C or N.")

    def test_alleles(self):
        # check for row duplications
        result = self.sqlContext.sql("SELECT chrom, pos, ref, alt, subject_id, sample_id FROM parquet_tbl WHERE \
            length(ref) = length(alt) AND alt <> '.' GROUP BY chrom, pos, ref, alt, subject_id, sample_id \
        HAVING COUNT(*) > 2")
        self.assertTrue(result.count() == 0, "More than two alleles found per locus.")

    def subjects_loaded(self):
        if str.lower(self.platform) == 'cgi':
            vcf_file = 'hdfs://nameservice1/vlad/wgs_cg.db/vcf_file_test/'
            file_tbl = self.sqlContext.read.parquet("{}".format(vcf_file))
            # register parquet file as temp table to query
            vcf_file_tbl = file_tbl.registerTempTable("vcf_file_tbl")
            result = self.sqlContext.sql("WITH a AS (SELECT DISTINCT subject_id FROM parquet_tbl) \
                                         SELECT a.*, b.subject_id \
                                         FROM a \
                                         FULL OUTER JOIN vcf_file_tbl b \
                                         ON a.subject_id = b.subject_id \
                                         WHERE (a.subject_id IS NULL or b.subject_id IS NULL)")
            self.assertTrue(result.count() == 0,
                            "All subject id's were not found in both tables.")

        elif str.lower(self.platform) == 'illumina':
            vcf_file = 'hdfs://nameservice1/vlad/wgs_ilmn.db/vcf_file_test/'
            file_tbl = self.sqlContext.read.parquet("{}".format(vcf_file))
            # register parquet file as temp table to query
            vcf_file_tbl = file_tbl.registerTempTable("vcf_file_tbl")
            result = self.sqlContext.sql("WITH a AS (SELECT DISTINCT subject_id FROM parquet_tbl) \
                                         SELECT a.*, b.subject_id \
                                         FROM a \
                                         FULL OUTER JOIN vcf_file_tbl b \
                                         ON a.subject_id = b.subject_id \
                                         WHERE (a.subject_id IS NULL or b.subject_id IS NULL)")
            self.assertTrue(result.count() == 0,
                        "All subject id's were not found in both tables.")

        else:
            print("Please enter either platform='illumina' or platform='cgi' in your class argument.")

    def chrom_per_subject(self):
        result = self.sqlContext.sql("with t as (select distinct subject_id, chrom \
                                     from vcf_variant) \
                                     select subject_id, count(*) as count from t \
                                     group by subject_id having count < 25")
        self.assertTrue(result.count() == 0,
                        "Not all chromosomes were loaded for the following subject_id's: \n {}".format(result))

    def tear_down(self):
        # clear memory cache
        self.sqlContext.clearCache()
        # close connection
        self.sc.stop()

    def run_tests(self):
        # check that table exists
        try:
            os.path.isfile(self.in_table)
        except Exception as e:
            print e

        # listing each test instead of using run_main
        # for ease of customization

        # check chrom, pos, sample_id column vals not empty
        self.check_empty("chrom")
        self.check_empty("ref")
        self.check_empty("subject_id")
        self.check_empty("sample_id")

        # check that all chroms loaded
        self.check_chrom_set()

        # check that each chrom contains at least 10000 rows
        self.check_chrom_count()

        # check that 'chr' prefix has been removed from chromosome
        self.check_chrom_prefix()

        # check that mitochondria has been renamed from MT to M
        self.check_chrom_mt()

         # check for no more than two alleles per person at each locus
        self.test_alleles()


###############
if __name__ == '__main__':

    # instantiate class
    spark = TestETL(in_table='hdfs://nameservice1/vlad/wgs_ilmn.db/vcf_variant/', is_vcf='yes',
                    platform='illumina', local_dir=os.getcwd(), hdfs_dir='/user/selasady/', appname='vcf_testing')

    if spark.is_vcf == 'yes':
        # run tests
        spark.run_tests()
        # check that all subjects listed in the vcf_file metadata were loaded
        spark.subjects_loaded()
        spark.chrom_per_subject()
        spark.tear_down()
    else:
        # run tests
        spark.run_tests()
        print ("Unit tests complete.")
        spark.tear_down()



# /itmi/wgs_cg.db/vcf_nocall
# /itmi/wgs_cg.db/vcf_variant