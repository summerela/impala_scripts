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
        self.conf = SparkConf().setAppName(self.appname) \
            .set("spark.sql.parquet.compression.codec", "snappy") \
            .set("spark.yarn.executor.memoryOverhead", 14336)
        #                    .set("spark.yarn.executor.cores", 1)
        self.sc = SparkContext(conf=self.conf)
        self.sqlC = SQLContext(self.sc)
        self.sqlC.sql("SET spark.sql.parquet.binaryAsString=true")
        self.sqlC.sql("SET spark.sql.parquet.cacheMetadata=true")
        #        self.sqlC.sql("SET spark.sql.parquet.compression.codec=snappy")
        self.in_table = "{}{}{}".format(self.spark_host_prefix, self.ilmn_db, 'vcf_distinct')
        #        self.var_df = self.sqlC.parquetFile(self.in_table)
        #        self.var_tbl = self.var_df.registerTempTable("var_tbl")
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
        files = '{}/chrom={}'.format(self.in_table, input_chrom)

        files = '/tmp/gv_test'
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
                      "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
                      "ANN[*].DISTANCE" "ANN[*].HGVS_C" \
                      "ANN[*].HGVS_P"'.format(self.snpsift_jar)

        rdd = self.sqlC.parquetFile(files) \
            .map(self.var2tsv) \
            .pipe(snpeff_cmd) \
            .pipe(self.snpeff_oneperline) \
            .pipe(extract_cmd) \
            .filter(lambda l: not l.startswith('#')).map(lambda l: Row(name=l))
            #  .map(lambda l: Row(chrom=tsv[0], pos=int(tsv[1], id=tsv[2])))



        df = self.sqlC.createDataFrame(rdd)
        print (df.take(5))
        # df.write.parquet('/tmp/gv_out')
        sys.exit(0)

    # select variants by chromosome
    #        var_df =  self.sqlC.sql("select * from var_tbl where chrom = '{}'".format(input_chrom))
    # run snpeff on query results
    #        snp_out = "{}/chr{}_snpeff.vcf".format(self.out_dir, input_chrom)
    #        snpeff_cmd = r'''java -d64 -Xmx32g -jar {snpeff} -t -v GRCh37.75 > {vcf_out}'''.format(snpeff=self.snpeff_jar, vcf_out=snp_out)

    # tsv_df = var_df.map(lambda r: r[0])
    # stuff= var_df.map(lambda line: line.split('\t')).collect()

    # var_df.rdd.map(lambda x: ",".join(map(str, x))).coalesce(1).saveAsTextFile("test.csv")



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