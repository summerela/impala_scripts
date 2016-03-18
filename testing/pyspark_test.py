#!/usr/bin/env python

from pyspark import SparkConf, SparkContext
conf = (SparkConf()
         .setMaster("glados15")
         .setAppName("testing spark")
         .set("spark.executor.memory", "1g"))
sc = SparkContext(conf = conf)

print (sc)






