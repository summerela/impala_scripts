#!/usr/bin/env pyspark

'''
Requirements:
    - spark 1.3 and above
    - python 2.7x
    - SnpEff 4.2 (build 2015-12-05)

Variables:
spark_host_prefix = hdfs prefix for cluster
in_files = folder containing variant parquet files
out_dir = directory to write snpeff parsed output

Purpose:
As part of the global variants pipeline, runs snpeff on variants located in the in_files directory

input = hdfs location of parquet formatted distinct variants created using gv_illumina.py in table called ilmn_vars
output = parquet formatted snpeff results parsed to one variant per line

'''

from pyspark import SparkContext, SparkConf, SQLContext
import os


class snpeff(object):

    # ITMI impala cluster
    tool_path = '/opt/cloudera/parcels/ITMI/'
    snpeff_jar = '{}share/snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}share/snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}share/snpEff/SnpSift.jar'.format(tool_path)

    def __init__(self, spark_host_prefix='hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/', in_files='/itmi/wgs_ilmn.df/ilmn_vars/',
                 out_dir = '/itmi/wgs_ilmn_new/snpeff_out'):
        self.spark_host_prefix = spark_host_prefix
        self.in_files = in_files
        self.out_dir = out_dir
        self.appname = "run_snpeff"
        self.conf = SparkConf().setAppName(self.appname) \
            .set("spark.sql.parquet.compression.codec", "snappy") \
        # disable yarn memory on itmi servers
            # .set("spark.yarn.executor.memoryOverhead", 14336)
        self.sc = SparkContext(conf=self.conf)
        self.sqlC = SQLContext(self.sc)
        self.sqlC.sql("SET spark.sql.parquet.binaryAsString=true")
        self.sqlC.sql("SET spark.sql.parquet.cacheMetadata=true")

    @staticmethod
    def check_outdir(output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    @staticmethod
    def var2tsv(r):
        '''
        Convert a variant RDD row into a tab delimited line in VCF format
        '''
        # RDD: chrom, pos,  var_id, ref, allele,       ...       stop, blk_pos
        # VCF: chrom, pos,  id,     ref, alt,    qual, filter, info
        return '\t'.join([str(r.chrom), str(r.pos), r.var_id, r.ref, r.allele,
                          '', '', ''])
    def run_snpeff(self):
        '''
        Run snpeff and parse output
        :return: vcf files of annotated variants for each chrom

        Test file is in /tmp/gv_test/chrom=1/blk_pos=1/ff*.parq.  Note that
        hadoop correctly interprets the chrom=1 and blk_pos=1 as a partitions
        and includes them as additional 'chrom' and blk_pos columns in the RDD

        '''


        # Do not run snpEff with stats or threads (Spark job is distributed and
        # threading is outside of resource management and likely won't even help)
        snpeff_cmd = 'java -d64 -Xmx16g -jar {} -noStats -v GRCh37.75'.format(self.snpeff_jar)

        extract_cmd = 'java -Xmx4g -jar {} extractFields - \
                      CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" \
                      "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].FEATURE" \
                      "ANN[*].FEATUREID" "ANN[*].BIOTYPE"\
                      "ANN[*].HGVS_C" "ANN[*].HGVS_P"'.format(self.snpsift_jar)

        rdd = self.sqlC.parquetFile(self.in_files) \
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
        df.write.parquet(self.out_dir)

    def shut_down(self):
        self.sqlC.clearCache()
        self.sc.stop()


############################################
if __name__ == "__main__":

    # instantiate instance (change depending on ITMI or ISB cluster)
    gv = snpeff(spark_host_prefix='localhost', in_files= '/itmi/wgs_ilmn.df/ilmn_vars/',
                out_dir='./snpeff_out')
    # create output directory if not exists
    gv.check_outdir(gv.out_dir)
    # run snpeff and parse output
    gv.run_snpeff()
    # shut down connection to spark
    gv.shut_down()