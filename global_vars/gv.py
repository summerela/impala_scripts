#!/usr/bin/env python

from pyspark import SparkContext, SparkConf, SQLContext
import subprocess as sp

class gv(object):

    # chroms = map(str, range(1, 23)) + ['X', 'Y', 'M']
    chroms = ['X', 'Y']
    var_blocks = range(0, 251)

    # ITMI impala cluster
    tool_path = '/opt/cloudera/parcels/ITMI/'
    snpeff_jar = '{}share/snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}share/snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}share/snpEff/SnpSift.jar'.format(tool_path)

    def __init__(self, spark_host_prefix='hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/', ilmn_db='wgs_ilmn.db/', \
                 anno_db='anno_grch37.db/'):
        self.spark_host_prefix = spark_host_prefix
        self.ilmn_db = ilmn_db
        self.anno_db = anno_db
        self.conf = SparkConf().setAppName("gv_pipeline")
        self.sc = SparkContext(conf=self.conf)
        self.sqlC = SQLContext(self.sc)
        self.sqlC.sql("SET spark.sql.parquet.binaryAsString=true")
        self.sqlC.sql("SET spark.sql.parquet.cacheMetadata=true")

    def register_table(self, prefix, in_table):
        print("Registering spark temp table {}...".format(in_table))
        in_df = self.sqlC.parquetFile("{}/{}".format(prefix, in_table))
        in_df.registerTempTable('{}'.format(in_table))

    @staticmethod
    def subprocess_cmd(command, cwd=os.getcwd()):
        '''
        Run programs in bash via subprocess
        :param command: command string as would be run on the command line
        :param input_dir: optional directory to run command in, default cwd
        :return: runs bash command
        '''
        print ("Running \n {}".format(command))
        ps = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd)
        try:
            print ps.communicate()
        except sp.CalledProcessError as e:
            print e

    def run_query(self, input_query):
        self.sqlC.sql(input_query)

    def shut_down(self):
        self.sqlC.clearCache()
        self.sc.stop()


############################################
if __name__ == 'main':

    gv = gv()


    gv.shut_down()
