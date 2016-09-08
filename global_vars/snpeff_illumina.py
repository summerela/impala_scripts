#!/usr/bin/env pyspark

from pyspark import SparkContext, SparkConf, SQLContext
from pyspark.sql import Row
import subprocess as sp
import os, csv, io, logging
import sys

logger = logging.getLogger('snpeff')
hdlr = logging.FileHandler('snpeff.log')
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)


class snpeff(object):
    # chroms = map(str, range(1, 23)) + ['X', 'Y', 'M']
    chroms = ['X']
    var_blocks = range(0, 251)

    # ITMI impala cluster
    tool_path = '/opt/cloudera/parcels/ITMI/'
    snpeff_jar = '{}share/snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}share/snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}share/snpEff/SnpSift.jar'.format(tool_path)

    def __init__(self, spark_host_prefix='hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/'):
        self.spark_host_prefix = spark_host_prefix
        self.appname = "run_snpeff"
        self.conf = SparkConf().setAppName(self.appname) \
            .set("spark.sql.parquet.compression.codec", "snappy") \
            .set("spark.yarn.executor.memoryOverhead", 14336)
        #                    .set("spark.yarn.executor.cores", 1)
        self.sc = SparkContext(conf=self.conf)
        self.sqlC = SQLContext(self.sc)
        self.sqlC.sql("SET spark.sql.parquet.binaryAsString=true")
        self.sqlC.sql("SET spark.sql.parquet.cacheMetadata=true")
        #        self.sqlC.sql("SET spark.sql.parquet.compression.codec=snappy")
        self.in_table = "{}/ilmn_db/ilmn_vars/".format(self.spark_host_prefix)


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

    def shut_down(self):
        self.sqlC.clearCache()
        self.sc.stop()

    @staticmethod
    def var2tsv(r):
        '''
        Convert a variant RDD row into a tab delimited line in VCF format
        '''
        # RDD: chrom, pos,  var_id, ref, allele,       ...       stop, blk_pos
        # VCF: chrom, pos,  id,     ref, alt,    qual, filter, info
        return '\t'.join([str(r.chrom), str(r.pos), r.var_id, r.ref, r.allele,
                          '', '', ''])

    def run_snpeff(self, input_chrom):
        '''
        Run snpeff by chromosome on ilmn_vars table
        :return: vcf files of annotated variants for each chrom
        '''
        # files = '{}/chrom={}'.format(self.in_table, input_chrom)
        files = '/itmi/wgs_ilmn_new/vcf_variant/chrom={}/'.format(input_chrom)

        # files = '/tmp/gv_test'
        # Test file is in /tmp/gv_test/chrom=1/blk_pos=1/ff*.parq.  Note that
        # hadoop correctly interprets the chrom=1 and blk_pos=1 as a partitions
        # and includes them as additional 'chrom' and blk_pos columns in the RDD

        # Do not run snpEff with stats or threads (Spark job is distributed and
        # threading is outside of resource management and likely won't even
        # help)
        snpeff_cmd = 'java -d64 -Xmx16g -jar {} -noStats -v GRCh37.75'.format(self.snpeff_jar)

        extract_cmd = 'java -Xmx4g -jar {} extractFields - \
                      CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" \
                      "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].FEATURE" \
                      "ANN[*].FEATUREID" "ANN[*].BIOTYPE"\
                      "ANN[*].HGVS_C" "ANN[*].HGVS_P"'.format(self.snpsift_jar)

        rdd = self.sqlC.parquetFile(files) \
            .map(self.var2tsv) \
            .pipe(snpeff_cmd) \
            .pipe(self.snpeff_oneperline) \
            .pipe(extract_cmd) \
            .filter(lambda l: not l.startswith('#')).map(lambda x: x.split('\t')) \
            .map(lambda l: (str(l[0]), int(l[1]), str(l[2]),
                            str(l[3]), str(l[4]), str(l[5]),
                            str(l[6]), str(l[7]), str(l[8]),
                            str(l[9]), str(l[10]), str(l[11]),
                            str(l[12]), str(l[13])
                            ))




        df = self.sqlC.createDataFrame(rdd)
        print (df.take(5))
        # df.write.parquet('/tmp/gv_out3')
        sys.exit(0)




############################################
if __name__ == "__main__":

    gv = snpeff(spark_host_prefix='localhost')

    gv.check_outdir('./snpeff_out')

    for chrom in gv.chroms:
        gv.run_snpeff(chrom)

gv.shut_down()