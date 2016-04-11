#!/usr/bin/env python

import os
import sys
import subprocess as sp
import findspark
findspark.init()

# # Path for spark source folder
# os.environ['SPARK_HOME']="/opt/cloudera/parcels/CDH/lib/spark"
# # Append pyspark  to Python Path
# sys.path.append("/opt/cloudera/parcels/CDH/lib/spark/python")

from pyspark import SparkContext, SparkConf

########################
### Connect to spark ###
########################

class spark(object):

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
        # for dependencies, add pyFiles=['file1.py','file2'] argument
        self.sc = SparkContext(conf=self.conf)

    def subprocess_cmd(self, command, local_dir):
        '''
        Run programs in bash via subprocess
        :param command: command string as would be run on the command line
        :param local_dir: optional directory to run command in, default cwd
        :return: runs bash command
        '''
        print ("Running \n {}".format(command))
        ps = sp.Popen(command, shell=True,stdout=sp.PIPE,stderr=sp.PIPE, cwd=local_dir)
        try:
           print ps.communicate()
        except sp.CalledProcessError as e:
             print e

    def hdfs_put(self, input_file):
        hdfs_mkdir_cmd = "hdfs dfs -mkdir {}".format(self.hdfs_dir)
        hdfs_put_cmd = "hdfs dfs -put {} {}".format(input_file, self.hdfs_dir)
        self.subprocess_cmd(hdfs_mkdir_cmd, os.getcwd())
        self.subprocess_cmd(hdfs_put_cmd, os.getcwd())

    def hdfs_read(self, input_dir):
        lines = self.sc.textFile(input_dir)
        return lines

###############
if __name__ == '__main__':

    spark_con = spark(os.getcwd(), '/user/selasady/testing3/', "local", "etl_test", 2)

    local_file = '/users/selasady/my_titan_itmi/impala_scripts/testing/test/tale_of_two_cities.txt'

    # write file to hdfs using sys call to hdfs put
    #spark_con.hdfs_put(local_file)
 
    # read in a file from hdfs
    print spark_con.hdfs_read(spark_con.hdfs_dir)



    # close connection
    spark_con.sc.stop()
    
