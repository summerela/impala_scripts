#!/usr/bin/env python

'''
Assumptions:
- python 2.7.10
- snpeff 4.2
- GRCh37.75
- impala ODBC driver accessed through python impyla module
- gv_illumina.sh script created table wgs_ilmn.illumina_vars
- set write permissions on output directory and hdfs directory
- cloudera impala cluster containing table of variants

Input:
impala table of variants to annotate

Output:
annotated variants uploaded to hdfs and converted
into an impala table

'''
import pandas as pd
from impala.dbapi import connect
import datetime
import subprocess as sp
from impala.util import as_pandas
import os

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

class snpeff_pipeline(object):

    # set related file paths

    # ITMI impala cluster
    tool_path = '/opt/cloudera/parcels/ITMI/'
    snpeff_jar = '{}share/snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}share/snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}share/snpEff/SnpSift.jar'.format(tool_path)

    # on ISB impala cluster
    # tool_path = '/users/selasady/my_titan_itmi/tools/'
    # snpeff_jar = '{}snpEff/snpEff.jar'.format(tool_path)
    # snpeff_oneperline = '{}snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    # snpsift_jar = '{}snpEff/SnpSift.jar'.format(tool_path)

    def __init__(self, out_dir, impala_host, impala_port, impala_user_name, hdfs_path):
        self.out_dir = out_dir
        self.impala_host = impala_host
        self.impala_port = impala_port
        self.impala_name = impala_user_name
        self.hdfs_path = hdfs_path
        self.chroms = map(str, range(1, 22)) + ['X', 'Y']
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

    def run_query(self, input_query):
        '''
        Query impala
        :param input_query: query to run as string
        :return: query results as pandas dataframe
        '''
        self.cur.execute(input_query)
        query_df = as_pandas(self.cur)
        return query_df

    def vars_to_snpeff(self):
        '''
        Run snpeff by chromosome on vcf_distinct table
        :return: vcf files of annoated variants for each chrom
        '''
        for chrom in self.chroms:
            # select variants by chromosome
            get_vars_query = "SELECT chrom as '#chrom', pos, var_id as id, ref, allele as alt, 100 as qual, \
                             'PASS' as filter, 'GT' as 'format', '.' as INFO from wgs_ilmn.illumina_vars \
                             where chrom = '{}'".format(chrom)
            var_df = self.run_query(get_vars_query)
            # run snpeff on query results
            if not var_df.empty:
                snp_out = "{}/chr{}_snpeff.vcf".format(self.out_dir, chrom)
                snpeff_cmd = r'''java -Xmx16g -jar {snpeff} -t GRCh37.75 > {vcf_out}'''.format(snpeff=self.snpeff_jar,
                                                                                                     vcf_out=snp_out)
                # run the subprocess command
                ps = sp.Popen(snpeff_cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=os.getcwd())
                try:
                    print ps.communicate(var_df.to_csv(sep='\t', header=True, index=False))
                except sp.CalledProcessError as e:
                    print e
            else:
                print("No variants found for chromosome {}".format(chrom))

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
        self.subprocess_cmd(parse_cmd, self.out_dir)

    def run_parse(self):
        '''
        parse snpeff output and create tsv files for
        upload to impala table
        :param out_dir: path to directory containing snpeff annotated vcf files
        :return: annotated tsv files for upload to impala
        '''
        for file in os.listdir(self.out_dir):
            if file.endswith('_snpeff.vcf'):
                self.parse_snpeff(file)

    ############################################
    ## Remove Header and add pos_block column ##
    ############################################

    def parse_tsv(self, input_tsv):
        '''
        remove header from tsv file and subset columns
        :param input_tsv: one-per line tsv file of snpeff annotated variants
        :return: final.tsv with no header and extraneous columns removed
        '''
        final_out = "{}/{}_final.tsv".format(self.out_dir, self.get_file_base(input_tsv, 1))
        final_df = pd.read_csv("{}/{}".format(self.out_dir, input_tsv), sep='\t', skiprows=1, header=None)
        # cant use seq in pandas df slicing
        final_df = final_df[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0]]
        final_df['pos_block'] = final_df[1].div(1000000).astype(int)
        final_df.to_csv(final_out, sep='\t', header=False, index=False)

    def run_parse_tsv(self):
        '''
        run the final tsv parsing function parse_tsv()
        :return: chrom_snpeff.tsv ready for upload to hdfs
        '''
        for file in os.listdir(self.out_dir):
            if file.endswith('_snpeff.tsv'):
                self.parse_tsv(file)

    #############################
    ## Upload results to hdfs  ##
    #############################

    def make_hdfs_dir(self):
        '''
        make a directory to upload snpeff tsv files
        :return: hdfs output directory
        '''
        mkdir_cmd = "hdfs dfs -mkdir {}".format(self.hdfs_out)
        self.subprocess_cmd(mkdir_cmd, self.out_dir)

    def upload_hdfs(self, in_tsv):
        '''
        upload chrom_snpeff.tsv files into hdfs_out
        :param in_tsv: chrom_snpeff.tsv files
        :return: files on hdfs directory
        '''
        upload_cmd = 'hdfs dfs -put {} {}'.format(in_tsv, self.hdfs_out)
        self.subprocess_cmd(upload_cmd, self.out_dir)

    def update_permissions(self):
        '''
        make sure we can write to hdfs directory
        :return: hdfs directory modified with read/write/view access to 777
        '''
        chown_cmd = "hdfs dfs -chmod 777 {}".format(self.hdfs_path)
        self.subprocess_cmd(chown_cmd, self.out_dir)

    def run_hdfs_upload(self):
        self.make_hdfs_dir()
        for file in os.listdir(self.out_dir):
            if file.endswith('_snpeff_final.tsv'):
                self.upload_hdfs(file)
        self.update_permissions()

    ##################
    ## Run routine  ##
    ##################

    def run_snpeff_routine(self):
        self.vars_to_snpeff()
        self.run_parse()
        self.run_parse_tsv()
        self.run_hdfs_upload()
        self.cur.close()


##########  Main Routine  ############
if __name__ == "__main__":

    # ITMI options
    out_dir = '/home/ec2-user/elasasu/impala_scripts/global_vars/gv_out'
    impala_host = 'glados14'
    impala_port = 21050
    impala_user_name = 'selasady'
    hdfs_path = '/user/selasady/'

    # ISB options
    # out_dir = '/titan/ITMI1/workspaces/users/selasady/impala_scripts/annotation/snpeff'
    # impala_host = 'glados14'
    # impala_port = 21050
    # impala_user_name = 'selasady'
    # hdfs_path = '/user/selasady/'

    #######################
    # run snpeff routines #
    #######################
    # instantiate class with user args
    snpeff = snpeff_pipeline(out_dir, impala_host, impala_port, impala_user_name, hdfs_path)
    # run the main routines
    snpeff.run_snpeff_routine()


