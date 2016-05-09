#!/usr/bin/env python

import pandas as pd
from impala.dbapi import connect
import time
import csv
import subprocess as sp
import numpy as np
from impala.util import as_pandas
import os

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

class run_snpeff(object):

    # set related file paths
    java = '/tools/java/jdk1.7/bin/java'
    gatk = '/users/selasady/my_titan_itmi/tools/GenomeAnalysisTK.jar'
    ref_fasta = '/users/selasady/my_titan_itmi/tools/human_g1k_v37.fasta'
    snpeff_jar = '/users/selasady/my_titan_itmi/tools/snpEff/snpEff.jar'
    snpeff_oneperline = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/vcfEffOnePerLine.pl'
    snpsift_jar = '/users/selasady/my_titan_itmi/tools/snpEff//SnpSift.jar'
    vcf_verify = '/users/selasady/my_titan_itmi/tools/snpEff/scripts/vcfBareBones.pl'

    def __init__(self, vcf_dir, impala_host, impala_port, impala_user_name):
        self.vcf_dir = vcf_dir
        self.impala_host = impala_host
        self.impala_port = impala_port
        self.impala_name = impala_user_name
        self.chroms = map(str, range(1, 22)) + ['X', 'Y', 'M']
        self.conn = connect(host=self.impala_host, port=self.impala_port, timeout=10000, user=self.impala_user)
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

    ######################################################
    # strip vcf formatting using snpeff vcfBareBones.pl ##
    ######################################################

    def filter_vcf(self, input_vcf):
        '''
        Use SnpEff VCF Bare Bones perl script to strip formatting from VCF file
        to minimize file size and compatibility issues
        :param input_vcf: path of vcf file to process
        :return: gzipped, trimmed vcf file file_trimmed.vcf.gz ready to feed to snpeff
        '''
        out_base = self.get_file_base(input_vcf, 2)
        strip_cmd = r'''{barebones} {} | gzip > {out_dir}/{out_vcf}_trimmed.vcf.gz'''.format(
            barebones=self.vcf_verify, out_dir=vcf_dir, out_vcf=out_base)
        self.subprocess_cmd(strip_cmd, self.vcf_dir)

    def run_filter(self, vcf_dir):
        for file in os.listdir(vcf_dir):
            if file.endswith('.vcf.gz'):
                self.subset_vcf(file)




##########  Main Routine  ############
vcf_dir = '/titan/ITMI1/workspaces/users/selasady/impala_scripts/annotation/snpeff'
impala_host = 'glados14'
impala_port = 21050
impala_user_name = 'selasady'

snpeff = run_snpeff(vcf_dir, impala_host, impala_port, impala_user_name)

snpeff.run_filter(vcf_dir)

