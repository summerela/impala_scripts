import os
import sys

# Path for spark source folder
os.environ['SPARK_HOME']="/opt/cloudera/parcels/CDH/lib/spark"

# Append pyspark  to Python Path
sys.path.append("/opt/cloudera/parcels/CDH/lib/spark/pythonx")

# try to import pyspark, throw error if..error
try:
    from pyspark import SparkContext
    from pyspark import SparkConf

    print ("Successfully imported Spark Modules")

except ImportError as e:
    print (e)