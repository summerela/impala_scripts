#!/usr/bin/env pyspark

import os
from pyspark import SparkContext, SparkConf, SQLContext
from pyspark.mllib.tree import RandomForest, RandomForestModel
from pyspark.mllib.util import MLUtils
from pyspark.sql.types import *

class pyspark_ml():

    def __init__(self, local_dir='./', hdfs_dir='/titan/ITMI1/workspaces/users/',
                 master='local', appname='spark_app', spark_mem=2):
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
        self.somatic_chr = sorted(map(int, range(1, 23)))
        self.sqlContext = SQLContext(self.sc)

    def tsv_to_libsvm(self, in_file):
        df = self.sqlContext.read.load(in_file,
                                  format='com.databricks.spark.csv',
                                  header='true',
                                  inferSchema='true')
        print df.describe().show()



    def random_forest(self, in_file):
        data = MLUtils.loadLibSVMFile(self.sc, in_file)
        # Split the data into training and test sets (30% held out for testing)
        (trainingData, testData) = data.randomSplit([0.7, 0.3])
        # Train a RandomForest model.
        #  Empty categoricalFeaturesInfo indicates all features are continuous.
        model = RandomForest.trainClassifier(trainingData, numClasses=2, categoricalFeaturesInfo={},
                                             numTrees=3, featureSubsetStrategy="auto",
                                             impurity='gini', maxDepth=4, maxBins=32)

        # Evaluate model on test instances and compute test error
        predictions = model.predict(testData.map(lambda x: x.features))
        labelsAndPredictions = testData.map(lambda lp: lp.label).zip(predictions)
        testErr = labelsAndPredictions.filter(lambda (v, p): v != p).count() / float(testData.count())
        print('Test Error = ' + str(testErr))
        print('Learned classification forest model:')
        print(model.toDebugString())


##########################
if __name__ == '__main__':

    # instantiate class with current wd, hdfs dir, spark evn, spark job name, num cores
    spark = pyspark_ml(os.getcwd(), '/user/selasady/', "local", "etl_test", 10)

    # specify name of input file on hdfs
    tsv_infile = 'sample_libsvm_data.txt'

    # parse vcf
    spark.tsv_to_libsvm(tsv_infile)

    # run tests
    # spark.random_forest(tsv_infile)