#!/usr/bin/env python

'''
Assumptions:
- python 2.7.10
- snpeff 4.2
- GRCh37.75
- impala accessible either directly on cluster or via impala ODBC/impyla
- gv_illumina.sh script created table wgs_ilmn.illumina_vars
- able to set write permissions on output directory and hdfs directory
- cloudera impala cluster containing table of variants
- Does not run on mitochondria

Input:
impala table of variants to annotate

Output:
annotated variants uploaded to hdfs and converted
into an impala table

'''
from subprocess import Popen, PIPE
import pandas as pd
from impala.dbapi import connect
import subprocess as sp
from impala.util import as_pandas
import os
import logging
import multiprocessing as mp


logger = logging.getLogger('snpeff')
hdlr = logging.FileHandler('snpeff.log')
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None


class snpeff(object):

    # chroms = map(str, range(2, 23)) + ['X', 'Y']
    chroms = ['X']
    var_blocks = range(0,251)

    # ITMI impala cluster
    tool_path = '/opt/cloudera/parcels/ITMI/'
    snpeff_jar = '{}share/snpEff/snpEff.jar'.format(tool_path)
    snpeff_oneperline = '{}share/snpEff/scripts/vcfEffOnePerLine.pl'.format(tool_path)
    snpsift_jar = '{}share/snpEff/SnpSift.jar'.format(tool_path)

    # ISB impala cluster
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
        self.conn = connect(host=self.impala_host, port=self.impala_port, timeout=10000, user=self.impala_name)
        self.cur = self.conn.cursor()
        self.hdfs_out = "{}/ilmn_snpeff".format(self.hdfs_path)

    # create function to run bash command with subprocess
    @staticmethod
    def subprocess_cmd(command, cwd=os.getcwd() ):
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

    ###################################################
    # Download Variants Table and run through snpeff ##
    ###################################################

    def pandas_query(self, input_query):
        '''
        Query impala and return results as pandas df
        :param input_query: query to run as string
        :return: query results as pandas dataframe
        '''
        try:
            print ("Running query: {}".format(input_query))
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


    def run_snpeff(self, input_chrom):
        '''
        Run snpeff by chromosome on ilmn_vars table
        :return: vcf files of annoated variants for each chrom
        '''
        # select variants by chromosome
        get_vars_query = "SELECT chrom as '#chrom', pos, var_id as id, ref, allele as alt, 100 as qual, \
                         'PASS' as filter, 'GT' as 'format', '.' as INFO from wgs_ilmn.ilmn_vars \
                         where chrom = '{}'".format(input_chrom)
        var_df = self.pandas_query(get_vars_query)
        # run snpeff on query results
        if not var_df.empty:
            snp_out = "{}/chr{}_snpeff.vcf".format(self.out_dir, input_chrom)
            snpeff_cmd = r'''java -d64 -Xmx4g -jar {snpeff} -t -v GRCh37.75 > {vcf_out}'''.format(snpeff=self.snpeff_jar,
                                                                                              vcf_out=snp_out)
            # run the subprocess command
            ps = sp.Popen(snpeff_cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=os.getcwd())
            stdout, stderr = ps.communicate(var_df.to_csv(sep='\t', header=True, index=False))
            if stdout:
                logger.info(stdout)
            if stderr:
                logger.error(stderr)
                print("Error encountered on chrom {}, check snpeff.log".format(input_chrom))

        else:
            print("No variants found for chromosome {}".format(input_chrom))

    ##########################################################
    ## Output SnpEff effects as tsv file, one effect per line ##
    ############################################################

# TODO remove snpeff file

    def parse_snpeff(self, input_chrom):
        in_name = "{}/chr{}_snpeff.vcf".format(self.out_dir, input_chrom)
        out_name = "{}/chr{}_snpeff.tsv".format(self.out_dir, input_chrom)
        parse_cmd = 'cat {vcf} | {perl} | java -Xmx4g -jar {snpsift} extractFields \
            - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].DISTANCE" \
            "ANN[*].HGVS_C" "ANN[*].HGVS_P" >> {out}'.format(vcf=in_name, perl=self.snpeff_oneperline, \
                                                            snpsift=self.snpsift_jar, out=out_name)
        self.subprocess_cmd(parse_cmd, self.out_dir)

    def parse_chunks(self, input_chrom):
        in_name = "{}/chr{}_snpeff.vcf".format(self.out_dir, input_chrom)
        chunksize = 100000
        reader = pd.read_csv(in_name, sep="\t", header=None, comment='#', chunksize=chunksize)
        pool = mp.Pool(4)
        for df in reader:
            pool.apply_async(self.parse_snpeff(df))

    def remove_splits(self, input_chrom):
        split_name = "{}/chr{}_snpeff.tsv".format(self.out_dir, input_chrom)
        if os.path.isfile(split_name):
            os.remove(split_name)
        else:
            raise SystemExit("The vcf file for chromosome {} was not split".format(input_chrom))






    ############################################
    ## Remove Header and add pos_block column ##
    ############################################

    def read_tsv_chunky(self, in_name):

        if os.path.isfile(in_name):
            reader = pd.read_csv("{}".format(in_name), sep='\t', comment='#', header=None,
                                   usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
                                   dtype={0: str, 1: str, 2: str, 3: str, 4: str, 5: str, 6: str, 7: str, 8: str,
                                          9: str, 10: str, 11: str, 12: str, 13: str, 14: int, 15: str},
                                   chunksize=100000)
            pool = mp.Pool(6)
            for chunk in reader:
                # reorder columns
                final_df = chunk[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0]]
                # add blk_pos column for processing
                final_df['blk_pos'] = final_df[1].div(1000000).astype(int)
                with open(final_out, 'a') as f:
                    final_df.to_csv(final_out, sep='\t', header=False, index=False)
                os.remove(in_name)
        else:
            raise SystemExit("A parsed snpeff tsv was not created for chromosome {}".format(in_name))


    def parse_tsv(self, input_chrom):
        '''
        remove header from tsv file and subset columns
        :param input_tsv: one-per line tsv file of snpeff annotated variants
        :return: final.tsv with no header and extraneous columns removed
        '''
        in_name = "{}/chr{}_snpeff.tsv".format(self.out_dir, input_chrom)
        final_out = "{}/chr{}_final.tsv".format(self.out_dir, input_chrom)




    #def parse_the_chunks(self, in_file):


    #############################
    ## Upload results to hdfs  ##
    #############################

    def make_hdfs_dir(self):
        '''
        make a directory to upload snpeff tsv files
        :return: hdfs output directory
        '''
        mkdir_cmd = "hdfs dfs -mkdir {}".format(self.hdfs_out)
        self.subprocess_cmd(mkdir_cmd)
        chown_cmd = "hdfs dfs -chmod 777 {}".format(self.hdfs_path)
        self.subprocess_cmd(chown_cmd)

    def upload_hdfs(self, input_chrom):
        '''
        upload chrom_snpeff.tsv files into hdfs_out
        :param in_tsv: chrom_snpeff.tsv files
        :return: files on hdfs directory
        '''
        up_file = "{}/chr{}_final.tsv".format(self.out_dir, input_chrom)
        upload_cmd = 'hdfs dfs -put {} {}'.format(up_file, self.hdfs_out)
        if os.path.isfile(up_file):
            self.subprocess_cmd(upload_cmd)
        else:
            raise SystemExit("A final tsv was not created for chromosome {}".format(input_chrom))

    def remove_final(self, input_chrom):
        '''
        check that final.tsv was uploaded to hdfs
        and delete from local dir
        :param input_chrom: chrom to check final.tsv on
        :return: removes final.tsv files from local dir
        '''
        remove_file = "{}/chr{}_final.tsv".format(self.out_dir, input_chrom)
        remove_cmd = "hdfs dfs -test -e {}".format(remove_file)
        process = Popen(remove_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        std_out, std_err = process.communicate()
        if std_err:
            raise SystemExit("Final tsv file was not uploaded for chromosome {}".format(input_chrom))
        else:
            os.remove(remove_file)

    ##################
    ## Run routine  ##
    ##################

    def run_snpeff_pipeline(self):
        # self.make_hdfs_dir()
        # self.check_outdir(self.out_dir)
        for chrom in snpeff.chroms:
            print ("Running snpeff on chromosome {} \n".format(chrom))
            self.run_snpeff(chrom)
            self.parse_chunks(chrom)
            # self.remove_splits(chrom)
            # self.parse_tsv(chrom)
            # self.upload_hdfs(chrom)
            # self.remove_final(chrom)
        self.cur.close()


# ##########  Main Routine  ############
if __name__ == "__main__":

    # ITMI options
    impala_host = 'localhost'
    impala_port = 21050
    impala_user_name = 'ec2-user'
    hdfs_path = 'elasasu/'
    out_dir = '/home/ec2-user/elasasu/impala_scripts/global_vars/illumina_gv'

    snp = snpeff(out_dir=out_dir, impala_host=impala_host, impala_port=impala_port, impala_user_name=impala_user_name,
                 hdfs_path=hdfs_path)
    snp.run_snpeff_pipeline()



########### unused snippets  #####################
# def split_snpeff(self, input_chrom):
#     '''
#     Parse snpeff output to contain one allele per row
#     :param input_vcf: snpeff annotation results
#     :return: tsv file for upload to impala table
#     '''
#     in_name = "{}/chr{}_snpeff.vcf".format(self.out_dir, input_chrom)
#     # split the snpeff file into smaller files for processing
#     split_cmd = "split -b 5G {}".format(in_name)
#     if os.path.isfile(in_name):
#         self.subprocess_cmd(split_cmd, self.out_dir)
#         os.remove(in_name)
#     else:
#         raise SystemExit("A vcf file was not created for chromosome {}".format(input_chrom))