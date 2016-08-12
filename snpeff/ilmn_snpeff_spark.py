#!/usr/bin/env pyspark

'''
Assumptions:
- python 2.7.10
- Spark 1.3.0 or abov
- snpeff 4.2
- GRCh37.75
- impala accessible either directly on cluster or via impala ODBC/impyla
- gv_illumina.sh script created table wgs_ilmn.illumina_vars
- able to set write permissions on output directory and hdfs directory
- cloudera impala cluster containing table of variants
- Does not run on mitochondria

Input:
impala table of variants to annotate

Output:
annotated variants uploaded to hdfs and converted
into an impala table

'''
from pyspark import SparkContext, SparkConf, SQLContext
from subprocess import Popen, PIPE
import pandas as pd
from impala.dbapi import connect
import subprocess as sp
from impala.util import as_pandas
import os
import logging
import multiprocessing as mp


logger = logging.getLogger('snpeff')
hdlr = logging.FileHandler('snpeff.log')
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None