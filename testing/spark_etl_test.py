#!/usr/bin/env pyspark

'''
Assumptions:
Python 2.6 or above
Spark 1.3.0 or above

Input:
Impala table in parquet format
Table contains at least these columns:
- chrom
- pos
- ref
- alt

If the vcf="yes" option is chosen, table must also contain these columns:
- subject_id
- sample_id

To Run:

- Edit main routine to instantiate class with correct table/options
- Type the following from command line: pyspark spark_etl_test.py  2>/dev/null

'''

import os
import unittest

from pyspark import SparkContext, SparkConf, SQLContext


############################
### Unit Test with Spark ###
############################

class TestETL(unittest.TestCase):

    def __init__(self, server_prefix='hdfs://nameservice1/itmi/', in_table = 'wgs_ilmn.db/vcf_variant', is_vcf='no', platform='illumina',
    local_dir='./', hdfs_dir='/itmi/', appname='etl_testing'):
        self.server_prefix = server_prefix
        self.in_table = in_table
        self.is_vcf = is_vcf
        self.platform = platform
        self.local_dir = local_dir
        self.hdfs_dir = hdfs_dir
        self.appname = appname
        self.conf = (SparkConf()
               .setAppName(self.appname))
        self.sc = SparkContext(conf=self.conf)
        self.chroms = sorted(map(str, range(1,23) + ["M", "X", "Y"]))
        self.sqlContext = SQLContext(self.sc)
        # encode any parquet binary formatting as strings
        self.sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")
        # cache table metadata for query speed
        self.sqlContext.sql("SET spark.sql.parquet.cacheMetadata=true")
        self.read_table = "%(server)s%(table)s" % {"server": self.server_prefix, "table": self.in_table}
        # self.data = self.sqlContext.parquetFile(self.read_table)
        self.data = self.sqlContext.read.parquet(self.read_table)
        # register parquet file as temp table to query
        self.parquet_tbl = self.data.registerTempTable("parquet_tbl")


    def check_empty(self, col_name):
        # Check that important columns do not contain null values
        result = self.sqlContext.sql("SELECT * FROM parquet_tbl WHERE '%(table)s' IS NULL" % {'table' : col_name}).count()
        self.assertTrue(result == 0, "Null values found in '%(column)s' column." % {'column' : col_name})

    def check_chrom_set(self):
        # check that all chromosomes were loaded
        result = self.sqlContext.sql("SELECT distinct chrom FROM parquet_tbl")
        uploaded_chroms = sorted(map(str, result.rdd.map(lambda r: r.chrom).collect()))
        self.assertEqual(self.chroms, uploaded_chroms,
                             "Not all chromosomes were uploaded.")

    def check_chrom_count(self):
        # check that sufficient data was uploaded for each chrom
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
        # Check that all subject id's were loaded by comparing table with metadata from redshift
        if str.lower(self.platform) == 'cgi':
            vcf_file = '%(server_prefix)s/wgs_cg.db/vcf_file/' % {'server_prefix': self.server_prefix}
            file_tbl = self.sqlContext.parquetFile("%(table)s" % {'table': vcf_file})
            # register parquet file as temp table to query
            vcf_file_tbl = file_tbl.registerTempTable("vcf_file_tbl")
            subjects = self.sqlContext.sql("SELECT DISTINCT subject_id FROM parquet_tbl").collect()
            result = self.sqlContext.sql("SELECT a.*, b.subject_id \
                                         FROM subjects a \
                                         FULL OUTER JOIN vcf_file_table b \
                                         ON a.subject_id = b.subject_id \
                                         WHERE (a.subject_id IS NULL or b.subject_id IS NULL)")
            self.assertTrue(result.count() == 0,
                            "All subject id's were not found in both tables.")

        elif str.lower(self.platform) == 'illumina':
            vcf_file = '%(server_prefix)s/wgs_ilmn.db/vcf_file/' % {'server_prefix': self.server_prefix}
            file_tbl = self.sqlContext.parquetFile("%(table)s" % {'table': vcf_file})
            # register parquet file as temp table to query
            vcf_file_tbl = file_tbl.registerTempTable("vcf_file_tbl")
            subjects = self.sqlContext.sql("SELECT DISTINCT subject_id FROM parquet_tbl").collect()
            result = self.sqlContext.sql("SELECT a.*, b.subject_id \
                                         FROM subjets a \
                                         FULL OUTER JOIN vcf_file_tbl b \
                                         ON a.subject_id = b.subject_id \
                                         WHERE (a.subject_id IS NULL or b.subject_id IS NULL)")
            self.assertTrue(result.count() == 0,
                        "All subject id's were not found in both tables.")

        else:
            print("Please enter either platform='illumina' or platform='cgi' in your class argument.")

    def chrom_per_subject(self):
        # make sure that all chromosomes were uploaded for each subject id
        subset_tbl = self.sqlContext.sql("SELECT DISTINCT subject_id, chrom from parquet_tbl").collect()
        result = self.sqlContext.sql("select subject_id, count(*) as count from subset_tbl \
                                     group by subject_id having count < 25")
        self.assertTrue(result.count() == 0,
                        "Not all chromosomes were loaded for the following subject_id's: \n '%(results)s" % {'results': result})

    def bracket_check(self):
        # check that structural variants and no calls were not loaded by searching for bracket notation
        result = self.sqlContext.sql("select * from parquet_tbl where instr(allele, '[') != 0 or instr(allele, ']') != 0")
        self.assertTrue(result.count() == 0, "Brackets were found in the allele field.")

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

        # list each test instead of using run_main
        # for ease of customization
        self.check_empty("chrom")
        self.check_empty("ref")
        self.check_chrom_set()
        self.check_chrom_count()
        self.check_chrom_prefix()
        self.check_chrom_mt()

###############
if __name__ == '__main__':

    # instantiate class
    spark = TestETL(server_prefix = 'hdfs://nameservice1/vlad/', in_table='wgs_ilmn.db/vcf_variant', is_vcf='yes',  \
    platform='illumina', local_dir=os.getcwd(), hdfs_dir='/user/selasady/', appname='etl_testing')

    if spark.is_vcf == 'yes':
        # run tests
        # spark.run_tests()
        # check that all subjects listed in the vcf_file metadata were loaded
        spark.subjects_loaded()
        spark.check_empty("subject_id")
        spark.check_empty("sample_id")
        spark.chrom_per_subject()
        spark.test_alleles()
        print ("Unit tests passed.")
        spark.tear_down()
    else:
        # run tests
        spark.run_tests()
        print ("Unit tests passed.")
        spark.tear_down()



# 'hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/wgs_ilmn.db/vcf_variant/'
# 'hdfs://nameservice1/vlad/wgs_ilmn.db/vcf_variant/'
# 'hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/wgs_cg.db/vcf_variant'


