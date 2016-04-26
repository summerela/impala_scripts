#!/usr/bin/env pyspark

import os
import unittest
import pandas as pd

from pyspark import SparkContext, SparkConf, SQLContext


############################
### Unit Test with Spark ###
############################

class etl_test(unittest.TestCase):

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
        self.somatic_chr = map( str, range(1,23) )
        self.sqlContext = SQLContext(self.sc)

        
    def tsv_to_df(self, input_file, header_bool, delim_var):
        # import file as dataframe, all cols will be imported as strings
        df = self.sqlContext.read.format("com.databricks.spark.csv").option("header", "true").option("delimiter", "\t").option("inferSchema", "true").load(input_file)
        # # cache df object to avoid rebuilding each time
        df.cache()
        return df

    def sparkDf_to_pandasDf(self, input_df):
        pandas_df = input_df.toPandas()
        return pandas_df

    def check_chrom_set(self, input_df):
        chrom_cols = input_df.select('chrom').distinct().collect()
        print ("The following chromosomes are loaded: {} \n").format(str(chrom_cols))
        assertItemsEqual(self.somatic_chr, chrom_cols, "Not all somatic chromosomes are loaded.")

    def check_chrom_empty(self, input_df):
        assert input_df.where(input_df.chrom.isNull()).count() == 0, "Null values found in chrom column."

    def check_ref_empty(self, input_df):
        assert input_df.where(input_df.ref.isNull()).count() == 0, "Null values found in ref column."

    def check_subject_empty(self, input_df):
        assert input_df.where(input_df.subject_id.isNull()).count() == 0, "Null values found in subject id column."

    def check_chrom_count(self, input_df):
        for name, group in pd_df['chrom'].groupby(pd_df['chrom']):
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

    def tear_down(self):
        # close connection
        self.sc.stop()



    #### tests for variant tables ###
    #


    # check that MT = M
    # if filter = pass, gt should never = '.'
    #



###############
if __name__ == '__main__':

    spark = etl_test(os.getcwd(), '/user/selasady/', "local", "etl_test", 2)

    file_df = spark.tsv_to_df('/user/selasady/etl_test/test_query.txt', 'True', '\t')
    pd_df = spark.sparkDf_to_pandasDf(file_df)

    ### chrom, pos, sample_id column vals not empty ###
    # spark.check_chrom_empty(file_df)
    # spark.check_ref_empty(file_df)
    # spark.check_subject_empty(file_df)

    ### check that all somatic chroms loaded ###
    # spark.check_chrom_set(file_df)

    ### check that each chrom contains at least 1000 rows
    # spark.check_chrom_count(pd_df)

    ### check that 'chr' prefix has been removed from chromosome
    spark.check_chrom_prefix(pd_df)

    ### ref cols contains expected values A T G C N and are uppercase ###
    # spark.ref_check(pd_df)








    
