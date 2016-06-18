#!/usr/bin/env pyspark

import os
from pyspark import SparkContext, SparkConf, SQLContext
from pyspark.mllib.linalg import Vectors
from pyspark.mllib.feature import StandardScaler
from pyspark.mllib.util import MLUtils
from pyspark.mllib.classification import SVMWithSGD
from pyspark.mllib.regression import LabeledPoint, LinearRegressionWithSGD
from pyspark.mllib.evaluation import RegressionMetrics
from pyspark.mllib.classification import LogisticRegressionWithLBFGS
from pyspark.mllib.classification import NaiveBayes
from pyspark.ml.feature import VectorAssembler
from pyspark.sql.functions import col




class pyspark_classify():

    def __init__(self, in_file, num_iterations=100, local_dir=os.getcwd(), hdfs_dir='/vlad/',
                 master='local', appname='spark_ml', spark_mem=4):
        self.in_file = in_file
        self.local_dir = local_dir
        self.hdfs_dir = hdfs_dir
        # number of times to iterate through training algorithm
        self.num_iterations = num_iterations
        # spark connection string options
        self.master = master
        self.appname = appname
        self.spark_mem = int(spark_mem)
        self.conf = (SparkConf()
                     .setMaster(self.master)
                     .setAppName(self.appname)
                     .set("spark.executor.memory", self.spark_mem))
        self.sc = SparkContext(conf=self.conf)
        # create sql context submodule for extra features
        self.sqlContext = SQLContext(self.sc)
        # read in parquet format data from hdfs
        self.data = self.sqlContext.read.parquet("{}".format(in_file))
        self.label = self.data.map(lambda x: x.label)
        self.features = self.data.map(lambda x: x[['qual', 'filter', 'rsid']])
        # create a set of labeled points to pass to svm model
        # self.labeled_points = self.data.map(lambda row: LabeledPoint(row[0], row[1:]))

    def scale_split(self):
        # create scaler model to transform data to mean of 0 and scale data by stdev
        scaler = StandardScaler(withMean=True, withStd=True).fit(self.features)
        # apply scaling to data features
        scalerModel = scaler.fit(self.data)
        scaledData = scalerModel.transform(self.data)
        # create training and test set
        trainingData, testData = scaledData.randomSplit([0.7, 0.3])
        # cach the training data for speed
        trainingData.cache()
        return trainingData, testData

    def getMSE(self, testingData, model):
        # get mean squared error for each model
        valuesAndPreds = testingData.map(lambda p: (p.label, model.predict(p.features[0])))
        MSE = valuesAndPreds.map(lambda (v, p): (v - p) ** 2).reduce(lambda x, y: x + y) / valuesAndPreds.count()
        return MSE

    def run_svm(self, training_data):
        # Build the model
        svm_model = SVMWithSGD.train(training_data, iterations=int(self.num_iterations))
        print("Mean Squared Error = " + str(self.getMSE(svm_model)))
        return svm_model

    def run_linear_regression(self, training_data):
        linear_regression_model = LinearRegressionWithSGD.train(training_data, iterations=self.num_iterations, step=0.1)
        return linear_regression_model

    def run_logistic_regression(self, training_data):
        logistic_regression_model = LogisticRegressionWithLBFGS.train(training_data, iterations=self.num_iterations)
        return logistic_regression_model

    def run_naive_bayes(self, training_data):
        bayes_model = NaiveBayes.train(training_data, lambda_=1.0)
        return bayes_model







##########################
if __name__ == '__main__':

    # instantiate class with input data file, num of iterations, current wd, hdfs dir, spark evn, spark job name, num cores
    spark = pyspark_classify('hdfs://nameservice1/vlad/wgs_ilmn.db/svm_parquet', 100, os.getcwd(), '/vlad/', "local", "etl_test", 10)


    # parsedData = spark.parse_data()
    #
    # print parsedData.take(2)






    # train, test = spark.scale_split()
    #
    # spark.run_svm(train)



    # def tsv_to_libsvm(self, in_file):
    #     df = self.sqlContext.read.load(in_file,
    #                               format='com.databricks.spark.csv',
    #                               header='false',
    #                               inferSchema='true',
    #                               comment = '#',
    #                               delimiter='\t' )
    #     return df

    # from pyspark.mllib.tree import RandomForest, RandomForestModel
    # def random_forest(self, in_file):
    #     data = MLUtils.loadLibSVMFile(self.sc, in_file)
    #     # Split the data into training and test sets (30% held out for testing)
    #     (trainingData, testData) = data.randomSplit([0.7, 0.3])
    #     # Train a RandomForest model.
    #     #  Empty categoricalFeaturesInfo indicates all features are continuous.
    #     model = RandomForest.trainClassifier(trainingData, numClasses=2, categoricalFeaturesInfo={},
    #                                          numTrees=3, featureSubsetStrategy="auto",
    #                                          impurity='gini', maxDepth=4, maxBins=32)
    #
    #     # Evaluate model on test instances and compute test error
    #     predictions = model.predict(testData.map(lambda x: x.features))
    #     labelsAndPredictions = testData.map(lambda lp: lp.label).zip(predictions)
    #     testErr = labelsAndPredictions.filter(lambda (v, p): v != p).count() / float(testData.count())
    #     print('Test Error = ' + str(testErr))
    #     print('Learned classification forest model:')
    #     print(model.toDebugString())