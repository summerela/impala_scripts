import os
import sys
import subprocess as sp


########################
### Connect to spark ###
########################
def subprocess_cmd(command, input_dir):
    '''
    Run programs in bash via subprocess
    :param command: command string as would be run on the command line
    :param input_dir: optional directory to run command in, default cwd
    :return: runs bash command
    '''
    print ("Running \n {}".format(command))
    ps = sp.Popen(command, shell=True,stdout=sp.PIPE,stderr=sp.PIPE, cwd=input_dir)
    try:
       print ps.communicate()
    except sp.CalledProcessError as e:
         print e

def hdfs_put(input_file, output_dir):
    hdfs_mkdir_cmd = "hdfs dfs -mkdir {}".format(output_dir)
    hdfs_put_cmd = "hdfs dfs -put {} {}".format(input_file, output_dir)
    subprocess_cmd(hdfs_mkdir_cmd, os.getcwd())
    subprocess_cmd(hdfs_put_cmd, os.getcwd())

# Path for spark source folder
os.environ['SPARK_HOME']="/opt/cloudera/parcels/CDH/lib/spark"

# Append pyspark  to Python Path
sys.path.append("/opt/cloudera/parcels/CDH/lib/spark/python")

# try to import pyspark, throw error if..error
from pyspark import SparkContext
from pyspark import SparkConf


conf = (SparkConf()
           .setMaster("local")
           .setAppName("test spark")
           .set("spark.executor.memory", "1g"))

# for dependencies, add pyFiles=['file1.py','file2'] argument
sc = SparkContext(conf=conf)

################
if __name__ == '__main__':

    # create new folder, wont write into existing
    hdfs_path = '/user/selasady/testing3/'
    local_file = '/users/selasady/my_titan_itmi/impala_scripts/testing/test/tale_of_two_cities.txt'

    # write file to hdfs using sys call to hdfs put
    hdfs_put(local_file, hdfs_path)   
 
    # testing pyspark by filtering out lines and counting non-empty
    lines = sc.textFile('/user/selasady/testing2/tale_of_two_cities.txt')
    non_empty_lines = lines.filter( lambda x: len(x) > 0)
    print (non_empty_lines.count())

    # word count map/reduce example
    words = non_empty_lines.flatMap(lambda x: x.split())
    wordcounts = words.map(lambda x: (x, 1)).reduceByKey(lambda x,y:x+y).map(lambda x:(x[1], x[0])).sortByKey(False)
    print (wordcounts.take(10))
    
    
