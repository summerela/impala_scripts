#!/usr/bin/env python

from impala.dbapi import connect
from impala.util import as_pandas
import subprocess as sp

def main():

    class connect_impala(object):

        # set related file paths

        def __init__(self, impala_host, impala_port, impala_user_name):
            self.impala_host = impala_host
            self.impala_port = impala_port
            self.impala_name = impala_user_name
            self.chroms = map(str, range(1, 22)) + ['X', 'Y']
            self.blk_pos = map(str, range(0, 250))
            self.conn = connect(host=self.impala_host, port=self.impala_port, timeout=10000, user=self.impala_name)
            self.cur = self.conn.cursor()

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
            try:
                self.cur.execute(input_query)
            except Exception as e:
                print e

################################
if __name__ == "__main__":
    main()