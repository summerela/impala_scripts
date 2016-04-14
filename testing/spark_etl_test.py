#!/usr/bin/env python

import os
import sys
import subprocess as sp


# Path for spark source folder
os.environ['SPARK_HOME']="/opt/cloudera/parcels/CDH/lib/spark"
# Append pyspark  to Python Path
sys.path.append("/opt/cloudera/parcels/CDH/lib/spark/python")

from pyspark import SparkContext, SparkConf

############################
### Unit Test with Spark ###
############################

class etl_test(unittest2.TestCase):

    def __init__(self, local_dir='./', hdfs_dir='/users/selasady/', master='master', appname='spark_job', spark_mem=2):
        self.local_dir = local_dir
        self.hdfs_dir = hdfs_dir
        self.master = master
        self.appname = appname
        self.spark_mem = int(spark_mem)
        self.conf = (SparkConf()
               .setMaster(self.master)
               .setAppName(self.appname)
               .set("spark.executor.memory", self.spark_mem))

    def set_up(self):
        # for dependencies, add pyFiles=['file1.py','file2'] argument
        self.sc = SparkContext(conf=self.conf)

    def tear_down(self):
        # close connection
        self.sc.stop()

    #### tests for variant tables ###
    # chrom set contains at least [1-22,X,Y,MT]
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

    spark = etl_test(os.getcwd(), '/user/selasady/testing3/', "local", "etl_test", 2)






    
