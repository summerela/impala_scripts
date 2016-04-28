#!/usr/bin/env pyspark

import os
import unittest
import urllib, json
import time

from pyspark import SparkContext, SparkConf, SQLContext


############################
### Unit Test with Spark ###
############################

class TestETL(unittest.TestCase):

    def __init__(self, local_dir='./', hdfs_dir='/titan/ITMI1/workspaces/users/', master='local', appname='spark_app', spark_mem=2):
        self.local_dir = local_dir
        self.hdfs_dir = hdfs_dir
        self.master = master
        self.appname = appname
        self.spark_mem = int(spark_mem)
        self.conf = (SparkConf()
               .setMaster(self.master)
               .setAppName(self.appname)
               .set("spark.executor.memory", self.spark_mem))
        self.sc = SparkContext(conf=self.conf)
        self.somatic_chr = sorted(map( int, range(1,23) ))
        self.sqlContext = SQLContext(self.sc)

        
    def tsv_to_df(self, input_file):
        # import file as dataframe, all cols will be imported as strings
        df = self.sqlContext.read.format("com.databricks.spark.csv").option("header", "true").option("delimiter", "\t").option("inferSchema", "true").load(input_file)
        # # cache df object to avoid rebuilding each time
        df.cache()
        # register as temp table for querying
        df.registerTempTable("spark_df")
        return df

    def sparkDf_to_pandasDf(self, input_df):
        pandas_df = input_df.toPandas()
        return pandas_df

    def check_chrom_empty(self, input_df):
        assert input_df.where(input_df.chrom.isNull()).count() == 0, "Null values found in chrom column."

    def check_ref_empty(self, input_df):
        assert input_df.where(input_df.ref.isNull()).count() == 0, "Null values found in ref column."

    def check_subject_empty(self, input_df):
        assert input_df.where(input_df.subject_id.isNull()).count() == 0, "Null values found in subject_id column."

    def check_allele_empty(self, input_df):
        assert input_df.where(input_df.allele.isNull()).count() == 0, "Null values found in allele column."

    def check_sample_empty(self, input_df):
        assert input_df.where(input_df.sample_id.isNull()).count() == 0, "Null values found in sample_id column."

    def check_chrom_set(self, input_df):
        # chrom_cols = sorted(input_df.select('chrom').distinct().collect()[1])
        chrom_cols = sorted([int(x) for x in input_df['chrom'].unique()])
        self.assertItemsEqual(self.somatic_chr, sorted(chrom_cols),
                              "Not all somatic chromosomes are loaded {}.".format(chrom_cols))

    def check_chrom_count(self, input_df):
        for name, group in input_df['chrom'].groupby(input_df['chrom']):
            self.assertGreater(group.count(), 999, "Chrom {} only has {} rows.".format(name, group.count()))

    def check_chrom_prefix(self, input_df):
        chroms = input_df['chrom'].unique()
        chrom_string = ''.join(str(e) for e in chroms)
        self.assertNotIn(chrom_string, 'chr')

    def check_chrom_mt(self, input_df):
        chroms = input_df['chrom'].unique()
        chrom_string = ''.join(str(e) for e in chroms)
        self.assertNotIn(chrom_string, 'mt')

    def ref_check(self, input_df):
        ref_col = input_df['ref'].unique()
        ref_set = set(''.join(sorted(ref_col)))
        result = [s for s in ref_set if s.strip('ATGCN')]
        assert len(result) == 0, "The following values were found in the ref field: {}".format(result)

    def test_alleles(self, input_df):
        no_indels = input_df[input_df['alt'].str.len() == input_df['ref'].str.len()]
        pos_group = no_indels.groupby(['subject_id', 'chrom', 'pos', 'ref'])
        problems = pos_group.filter(lambda x: len(x) > 2)
        assert len(problems) == 0, "Found more than two alleles for the following subject/position: \n {}".format(problems)

    def test_ref(self, input_df):
        ref_set = input_df[['chrom', 'pos', 'ref']][(input_df['ref'] != '.') and (input_df[input_df['alt'].str.len() == input_df['ref'].str.len()])]
        ref_group = ref_set.groupby(['chrom', 'pos'])
        for name, group in ref_group:
            assert len(set(group['ref'])) < 2, "Different reference alleles were found at {}".format(name)

    def get_ref_snp(self, input_df):
        random_snp = input_df[['chrom', 'pos', 'ref', 'alt']][input_df['ref'] != '.'].sample(1)
        return random_snp

    def query_snp(self, input_snp):
        exclude = ['ins', 'del']
        query_string = "chr{}:{}".format(input_snp.iloc[0]['chrom'], input_snp.iloc[0]['pos'])
        test_url = 'http://myvariant.info/v1/query?q={}'.format(query_string)
        response = urllib.urlopen(test_url)
        data = json.loads(response.read())
        ref_var = []
        while len(ref_var) == 0:
            for x in range(0,4):
                try:
                     ref_hit = data['hits'][0]['_id'] if len(data['hits'][0]['_id']) > 0 else 'null'
                     if any(x in exclude for x in ref_hit):
                         pass
                except:
                     time.sleep(10)
            else:
                ref_var.append(ref_hit)
                ref_out = ''.join(ref_var)
                return ref_out

    def check_one_based(self, input_df):
        test_snp = self.get_ref_snp(input_df)
        ref_snp = self.query_snp(test_snp)
        file_snp =  "chr{}:g.{}{}>{}".format(test_snp.iloc[0]['chrom'], test_snp.iloc[0]['pos'],  test_snp.iloc[0]['ref'],  test_snp.iloc[0]['alt'])
        self.assertEqual(ref_snp, file_snp, "File variant {} does not match reference {}".format(file_snp, ref_snp))

    def tear_down(self):
        # close connection
        self.sc.stop()

    def run_tests(self, input_df):

        try:
            os.path.isfile(input_df)
        except Exception as e:
            print e

        file_df = self.tsv_to_df(input_df)

        pd_df = self.sparkDf_to_pandasDf(file_df)

        # check chrom, pos, sample_id column vals not empty
        self.check_chrom_empty(file_df)
        self.check_ref_empty(file_df)
        self.check_subject_empty(file_df)
        self.check_allele_empty(file_df)
        self.check_sample_empty(file_df)

        # check that all somatic chroms loaded
        self.check_chrom_set(pd_df)

        # check that each chrom contains at least 1000 rows
        self.check_chrom_count(pd_df)

        # check that 'chr' prefix has been removed from chromosome
        self.check_chrom_prefix(pd_df)

        # ref cols contains expected values A T G C N and are uppercase ###
        self.ref_check(pd_df)

        # check for no more than two alleles per person at chrom, pos, ref
        self.test_alleles(pd_df)

        # # ref allele should be the same at each position
        self.test_ref(pd_df)

        # verify 1-based coords
        # self.check_one_based(pd_df)

        self.tear_down()





###############
if __name__ == '__main__':

    # instantiate class with current wd, hdfs dir, spark evn, spark job name, num cores
    spark = TestETL(os.getcwd(), '/user/selasady/', "local", "etl_test", 10)


    # specify input file to process
    tsv_infile = '/user/selasady/test_set.txt'

    # run tests
    spark.run_tests(tsv_infile)






    
