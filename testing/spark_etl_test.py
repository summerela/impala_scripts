#!/usr/bin/env pyspark

import os
import unittest
import pandas as pd

from pyspark import SparkContext, SparkConf, SQLContext
from pyspark.sql import Row

############################
### Unit Test with Spark ###
############################

class etl_test():

    def __init__(self, local_dir='./', hdfs_dir='/users/selasady/', master='local', appname='spark_app', spark_mem=2):
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
        lines = self.sc.textFile(input_file)
        parts = lines.map(lambda line: line.split('\t'))
        header = lines.first().split()
        header_map = ','.join(["p[{}]".format(header.index(s)) for s in header])
        stuff = parts.map(lambda p :(header_map).strip())
        # fields = [self.sqlContext.StructField(field_name, StringType(), True) for field_name in header.split()]
        # schema = self.sqlContext.StructField(fields)
        print stuff
        # schema = ','.join(["{} = line[ {} ]".format(s, header.index(s)) for s in header])


    def check_chrom(self, input_file):
        input_df = self.tsv_to_df(input_file)
        rdd_df = self.sqlContext.read(input_df)
        rdd_df.printSchema()


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

    spark.tsv_to_df('/user/selasady/etl_test/test_query.txt')






    
