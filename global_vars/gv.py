#!/usr/bin/env pyspark

from pyspark import SparkContext, SparkConf, SQLContext
import subprocess as sp
import os, csv, io, logging
import pandas as pd

logger = logging.getLogger('snpeff')
hdlr = logging.FileHandler('snpeff.log')
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

class snpeff(object):

    # chroms = map(str, range(1, 23)) + ['X', 'Y', 'M']
    chroms = ['X', 'Y']
    var_blocks = range(0, 251)

    # ITMI impala cluster
    tool_path = '/opt/cloudera/parcels/ITMI/'
    snpeff_jar = '{}share/snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}share/snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}share/snpEff/SnpSift.jar'.format(tool_path)

    def __init__(self, spark_host_prefix='hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/', ilmn_db='wgs_ilmn.db/', \
                 anno_db='anno_grch37.db/', out_dir='./'):
        self.spark_host_prefix = spark_host_prefix
        self.ilmn_db = ilmn_db
        self.anno_db = anno_db
        self.appname = "run_snpeff"
        self.conf = SparkConf().setAppName(self.appname)
        self.sc = SparkContext(conf=self.conf)
        self.sqlC = SQLContext(self.sc)
        self.sqlC.sql("SET spark.sql.parquet.binaryAsString=true")
        self.sqlC.sql("SET spark.sql.parquet.cacheMetadata=true")
        self.in_table = "{}{}{}".format(self.spark_host_prefix, self.ilmn_db, 'vcf_distinct')
        self.var_df = self.sqlC.parquetFile(self.in_table)
        self.var_tbl = self.var_df.registerTempTable("var_tbl")
        self.out_dir = out_dir

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

    @staticmethod
    def check_outdir(output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def run_query(self, input_query):
        self.sqlC.sql(input_query)

    def shut_down(self):
        self.sqlC.clearCache()
        self.sc.stop()

    def run_snpeff(self, input_chrom):
        '''
        Run snpeff by chromosome on ilmn_vars table
        :return: vcf files of annoated variants for each chrom
        '''
        # select variants by chromosome
        var_df =  self.sqlC.sql("select * from var_tbl where chrom = '{}'".format(input_chrom))
        # run snpeff on query results
        snp_out = "{}/chr{}_snpeff.vcf".format(self.out_dir, input_chrom)
        snpeff_cmd = r'''java -d64 -Xmx32g -jar {snpeff} -t -v GRCh37.75 > {vcf_out}'''.format(snpeff=self.snpeff_jar,
                                                                                          vcf_out=snp_out)

        stuff= var_df.map(lambda line: line.split('\t')).collect()
        print stuff.take(5)

        #stuff
        # run the subprocess command
        # ps = sp.Popen(snpeff_cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=os.getcwd())
        # stdout, stderr = ps.communicate(var_df)
        # if stdout:
        #     logger.info(stdout)
        # elif stderr:
        #     logger.error(stderr)
        #     raise SystemExit("Error encountered on chrom {}, check snpeff.log".format(input_chrom))



############################################
if __name__ == "__main__":

    gv = snpeff()

    gv.check_outdir('./snpeff_out')

    for chrom in gv.chroms:
        gv.run_snpeff(chrom)

    gv.shut_down()
