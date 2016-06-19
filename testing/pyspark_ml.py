#!/usr/bin/env pyspark

'''
 pyspark --packages com.databricks:spark-csv_2.10:1.2.0 pyspark_ml.py  2>/dev/null
'''

import os
from pyspark import SparkContext, SparkConf, SQLContext
from pyspark.mllib.feature import StandardScaler
from pyspark.mllib.classification import SVMWithSGD, LogisticRegressionWithLBFGS, NaiveBayes
from pyspark.mllib.regression import LabeledPoint, LinearRegressionWithSGD



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
        # read in data from hdfs
        self.data = self.sc.textFile("{}".format(self.in_file)).map(lambda line:[float(x) for x in line.split(' ')])
        self.features = self.data.map(lambda row: row[1:])
        self.labels = self.data.map(lambda row: row[0])

    def scale_data(self):
        # scale data to have 0 mean and std
        scaler = StandardScaler(withStd=True, withMean=False)
        scalerModel = scaler.fit(self.features)
        # Normalize each feature to have unit standard deviation.
        scaledFeatures = scalerModel.transform(self.features)
        transformedData = self.labels.zip(scaledFeatures)
        transformedData = transformedData.map(lambda row: LabeledPoint(row[0], [row[1]]))
        return transformedData

    def split_data(self):
        # split the data into training and testing sets
        scaled_data = self.scale_data()
        train_set, test_set = scaled_data.randomSplit([0.8,0.2])
        return train_set, test_set

    def getMSE(self, model, test_set):
        # use testing data to get mean squared error for each model
        valuesAndPreds = test_set.map(lambda p: (p.label, model.predict(p.features[0])))
        MSE = valuesAndPreds.map(lambda (v, p): (v - p) ** 2).reduce(lambda x, y: x + y) / valuesAndPreds.count()
        return MSE

    def getTrainingError(self, model, test_set):
        labelsAndPreds = test_set.map(lambda p: (p.label, model.predict(p.features)))
        trainErr = labelsAndPreds.filter(lambda (v, p): v != p).count() / float(test_set.count())
        return trainErr

    def run_svm(self, training_data):
        # Build the model
        svm_model = SVMWithSGD.train(training_data, iterations=int(self.num_iterations))
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

    def train_models(self):
        train, test = spark.split_data()
        svm = self.run_svm(train)
        svm_test = self.getTrainingError(svm, test)
        linear = self.run_linear_regression(train)
        linear_test = self.getMSE(linear, test)
        logistic = self.run_logistic_regression(train)
        logistic_test = self.getTrainingError(logistic, test)
        bayes = self.run_naive_bayes(train)
        bayes_test = self.getMSE(bayes, test)
        return svm_test, linear_test, logistic_test, bayes_test


##########################
if __name__ == '__main__':

    # instantiate class with input data file, num of iterations, current wd, hdfs dir, spark evn, spark job name, num cores
    # spark = pyspark_classify('hdfs://nameservice1/vlad/wgs_ilmn.db/svm_parquet', 100, os.getcwd(), '/vlad/', "local", "etl_test", 10)
    spark = pyspark_classify('sample_svm_data.txt', 100, os.getcwd(), '/User/selasady/', "local",
                             "etl_test", 10)

    svm_test, linear_test, logistic_test, bayes_test = spark.train_models()

    print "SVM training error = {}".format(svm_test)
    print "Linear regression mean squared error = {}".format(linear_test)
    print "Logistic regression training error = {}".format(logistic_test)
    print "Naive bayes mean standard error = {}".format(bayes_test)




















    ####################################################################################################################

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


    # self.data = self.sqlContext.read.parquet("{}".format(in_file))
    # self.label = self.data.map(lambda x: x.label)
    # self.features = self.data.map(lambda x: x[['qual', 'filter', 'rsid']])
    # create a set of labeled points to pass to svm model
    # self.labeled_points = self.data.map(lambda row: LabeledPoint(row[0], row[1:]))

    # create scaler model to transform data to mean of 0 and scale data by stdev
    # scaler = StandardScaler(withMean=True, withStd=True)
    # apply scaling to data features
    # model = scaler.fit(self.data)
    # result = model.transform(self.data)
    # create training and test set
    # trainingData, testData = scaledData.randomSplit([0.7, 0.3])
    # # cach the training data for speed
    # trainingData.cache()
    # return trainingData, testData
    # return result