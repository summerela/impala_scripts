#!/usr/bin/env pyspark

import os
import unittest
import pandas as pd

from pyspark import SparkContext, SparkConf, SQLContext,
import pyspark.sql.functions import fxn


############################
### Unit Test with Spark ###
############################

class etl_test():

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

        
    def tsv_to_df(self, input_file):
        # import file as dataframe, all cols will be imported as strings
        df = self.sqlContext.read.format("com.databricks.spark.csv").option("header", "true").option("delimiter", "\t").option("inferSchema", "true").load(input_file)
        # cache df object to avoid rebuilding each time
        df.cache()
        return df


    def check_chrom(self, input_df):
        chrom_cols = input_df.agg(self.sqlContext.(input_df.chrom)).collect()
        print chrom_cols

    def tear_down(self):
        # close connection
        self.sc.stop()



    #### tests for variant tables ###
    # chrom set contains at least [1-22]
    # chrom column vals not empty
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

    spark = etl_test(os.getcwd(), '/user/selasady/etl_test/', "local", "etl_test", 2)

    file_df = spark.tsv_to_df('/user/selasady/etl_test/test_query.txt')
    spark.check_chrom(file_df)






    
