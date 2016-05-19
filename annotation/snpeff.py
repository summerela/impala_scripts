#!/usr/bin/env python

import pandas as pd
from impala.dbapi import connect
import datetime
import subprocess as sp
from impala.util import as_pandas
import os

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

class run_snpeff(object):

    # set related file paths
    tool_path = '/users/selasady/my_titan_itmi/tools/'
    gatk = '{}GenomeAnalysisTK.jar'.format(tool_path)
    ref_fasta = '{}human_g1k_v37.fasta'.format(tool_path)
    snpeff_jar = '{}snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}snpEff/SnpSift.jar'.format(tool_path)
    vcf_verify = '{}snpEff/scripts/vcfBareBones.pl'.format(tool_path)

    def __init__(self, vcf_dir, impala_host, impala_port, impala_user_name, hdfs_path):
        self.vcf_dir = vcf_dir
        self.impala_host = impala_host
        self.impala_port = impala_port
        self.impala_name = impala_user_name
        self.hdfs_path = hdfs_path
        self.chroms = map(str, range(1, 22)) + ['X', 'Y']
        self.blk_pos = map(str, range(0, 250))
        self.conn = connect(host=self.impala_host, port=self.impala_port, timeout=10000, user=self.impala_name)
        self.cur = self.conn.cursor()
        self.now = datetime.datetime.now()
        self.today = str(self.now.strftime("%Y%m%d"))
        self.hdfs_out = "{}/snpeff_{}".format(self.hdfs_path, self.today)

    @staticmethod
    # create function to run bash command with subprocess
    def subprocess_cmd(command, input_dir):
        '''
        Run programs in bash via subprocess
        :param command: command string as would be run on the command line
        :param input_dir: optional directory to run command in, default cwd
        :return: runs bash command
        '''
        print ("Running \n {}".format(command))
        ps = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, cwd=input_dir)
        try:
            print ps.communicate()
        except sp.CalledProcessError as e:
            print e

    @staticmethod
    def get_file_base(in_file, num_to_remove):
        basename = str('.'.join(in_file.split('.')[:-int(num_to_remove)]) if '.' in in_file else in_file)
        return basename

    ###################################################
    # Download Variants Table and run through snpeff ##
    ###################################################

    def pandas_query(self, input_query):
        try:
            self.cur.execute(input_query)
            query_df = as_pandas(self.cur)
            return query_df
        except Exception as e:
            print e

    def run_query(self, input_query):
        print ("Running query: {}").format(input_query)
        try:
            self.cur.execute(input_query)
        except Exception as e:
            print e

    def vars_to_snpeff(self):
        for chrom in self.chroms:
            get_vars_query = "SELECT chrom as '#chrom', pos, var_id as id, ref, allele as alt, 100 as qual, \
                             'PASS' as filter, 'GT' as 'format', '.' as INFO from wgs_ilmn.ilmn_vars \
                             where chrom = '{}'".format(chrom)
            var_df = self.pandas_query(get_vars_query)
            if not var_df.empty:
                snp_out = "{}/chr{}_snpeff.vcf".format(self.vcf_dir, chrom)
                snpeff_cmd = r'''java -Xmx16g -jar {snpeff} -t GRCh37.75 > {vcf_out}'''.format(snpeff=self.snpeff_jar,
                                                                                                     vcf_out=snp_out)
                # run the subprocess command
                ps = sp.Popen(snpeff_cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=os.getcwd())
                try:
                    print ps.communicate(var_df.to_csv(sep='\t', header=True, index=False))
                except sp.CalledProcessError as e:
                    print e

    ##########################################################
    ## Output SnpEff effects as tsv file, one effect per line ##
    ############################################################
    def parse_snpeff(self, input_vcf):
        '''
        Parse snpeff output to contain one allele per row
        :param input_vcf: snpeff annotation results
        :return: tsv file for upload to impala table
        '''
        out_name = "{}.tsv".format(self.get_file_base(input_vcf, 1))
        parse_cmd = 'cat {vcf} | {perl} | java -Xmx16g -jar {snpsift} extractFields \
            - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].DISTANCE" \
            "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {out}'.format(vcf=input_vcf, perl=self.snpeff_oneperline, \
                                                         snpsift=self.snpsift_jar, out=out_name)
        self.subprocess_cmd(parse_cmd, self.vcf_dir)

    def run_parse(self):
        '''
        parse snpeff output and create tsv files for
        upload to impala table
        :param vcf_dir: path to directory containing snpeff annotated vcf files
        :return: annotated tsv files for upload to impala
        '''
        for file in os.listdir(self.vcf_dir):
            if file.endswith('_snpeff.vcf'):
                self.parse_snpeff(file)

    ############################################
    ## Remove Header and add pos_block column ##
    ############################################

    def parse_tsv(self, input_tsv):
        final_out = "{}/{}_final.tsv".format(self.vcf_dir, self.get_file_base(input_tsv, 1))
        final_df = pd.read_csv("{}/{}".format(self.vcf_dir, input_tsv), sep='\t', skiprows=1, header=None)
        # cant use seq in pandas df slicing
        final_df = final_df[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0]]
        final_df['pos_block'] = final_df[1].div(1000000).astype(int)
        final_df.to_csv(final_out, sep='\t', header=False, index=False)

    def run_parse_tsv(self):
        for file in os.listdir(self.vcf_dir):
            if file.endswith('_snpeff.tsv'):
                self.parse_tsv(file)

    #############################
    ## Upload results to hdfs  ##
    #############################

    def make_hdfs_dir(self):
        mkdir_cmd = "hdfs dfs -mkdir {}".format(self.hdfs_out)
        self.subprocess_cmd(mkdir_cmd, self.vcf_dir)

    def upload_hdfs(self, in_tsv):
        upload_cmd = 'hdfs dfs -put {} {}'.format(in_tsv, self.hdfs_out)
        self.subprocess_cmd(upload_cmd, self.vcf_dir)

    def update_permissions(self):
        chown_cmd = "hdfs dfs -chown -R impala:supergroup {}".format(self.hdfs_path)
        self.subprocess_cmd(chown_cmd, self.vcf_dir)

    def run_hdfs_upload(self):
        self.make_hdfs_dir()
        for file in os.listdir(self.vcf_dir):
            if file.endswith('_snpeff_final.tsv'):
                self.upload_hdfs(file)
        self.update_permissions()

    ###################
    ## Main routine  ##
    ###################
    def run_pipeline(self):
        print ("Running SnpEff on variants... \n")
        self.vars_to_snpeff()
        print ("Parsing snpeff output... \n")
        self.run_parse()
        self.run_parse_tsv()
        print ("Uploading snepff files to HDFS")
        self.run_hdfs_upload()


##########  Main Routine  ############
if __name__ == "__main__":
    run_snpeff.run_pipeline()


