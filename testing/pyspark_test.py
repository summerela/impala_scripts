#!/usr/bin/env python

from pyspark import SparkConf, SparkContext

conf = (SparkConf()
     .setMaster("glados14")
     .setAppName("testing spark")
     .set("spark.executor.memory", "1g")
     .setExecutorEnv('PYTHONPATH')

sc = SparkContext(conf=conf)