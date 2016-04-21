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


    def check_chrom_set(self, input_df):
        chrom_cols = input_df.select('chrom').distinct().collect()
        print ("The following chromosomes are loaded: {} \n").format(str(chrom_cols))
        self.assertItemsEqual(self.somatic_chr, chrom_cols, "Not all somatic chromosomes are loaded.")

    def check_chrom_empty(self, input_df):
        assert input_df.where(input_df.chrom.isNull()).count() == 0, "Null values found in chrom column."

    def check_chrom_count(self, input_df):
        print input_df.groupBy('chrom').count().show()

    def tear_down(self):
        # close connection
        self.sc.stop()



    #### tests for variant tables ###
    #
    # chrom, pos, sample_id column vals not empty
    # chrom cols contain at least one A T G C
    # gt contains at least one of 0/1, 1/1, 1/2
    # if gt = 0/0 at this chrom/pos allele_idx = (etc) ?
    # zygosity column is not all null
    # check that chr prefix has been removed
    # check that MT = M
    # if filter = pass, gt should never = '.'
    #



###############
if __name__ == '__main__':

    spark = etl_test(os.getcwd(), '/user/selasady/', "local", "etl_test", 2)

    file_df = spark.tsv_to_df('/user/selasady/etl_test/test_query.txt', 'True', '\t')

    # spark.check_chrom_set(file_df)
    # spark.check_chrom_empty(file_df)
    spark.check_chrom_count(file_df)







    
