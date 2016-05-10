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

class snpeff_pipeline(object):

    # set related file paths
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
        strip_cmd = r'''zcat {vcf} | {barebones} | gzip > {out_dir}/{out_vcf}_trimmed.vcf.gz'''.format(
            vcf=input_vcf, barebones=self.vcf_verify, out_dir=vcf_dir, out_vcf=out_base)
        self.subprocess_cmd(strip_cmd, self.vcf_dir)

    def run_filter(self, vcf_dir):
        '''
        run filter_vcf() on input directory
        :param vcf_dir: path to directory of vcf files to parse
        :return: _trimmed.vcf.gz files to run through snpeff
        '''
        for file in os.listdir(vcf_dir):
            if file.endswith('.vcf.gz'):
                self.filter_vcf(file)

    ###########################################
    # run _trimmed vcf files through snpeff  ##
    ###########################################

    def snpeff(self, input_vcf):
        '''
        Runs _trimmed.vcf.gz files created with filter_vcf() function
        through snpeff with the following options:
            Xmx16g= use 16g of memory for java
            t = Use multiple threads
        :param input_vcf: _trimmed.vcf.gz vcf file created with filter_vcf()
        :return: snpeff.vcf file with annotated snps
        '''
        snp_out = "{}_snpeff.vcf".format(self.get_file_base(input_vcf, 2))
        snpeff_cmd = r'''java -Xmx16g -jar {snpeff} closest -t GRCh37.75 {vcf} > {vcf_out}'''.format(
            snpeff = self.snpeff_jar, vcf= input_vcf, vcf_out= snp_out)
        self.subprocess_cmd(snpeff_cmd, self.vcf_dir)

    def run_snpeff(self, vcf_dir):
        '''
        run snpeff on all _trimmed.vcf.gz files in vcf_dir
        :param vcf_dir: path to _trimmed.vcf.gz files to run through snpeff
        :return: snpeff annotated vcf files
        '''
        for file in os.listdir(vcf_dir):
            if file.endswith('_trimmed.vcf.gz'):
                self.snpeff(file)

    ##########################################################
    ## Output SnpEff effects as tsv file, one effect per line ##
    ############################################################
    def parse_snpeff(self, input_vcf):
        '''
        Parse snpeff output to contain one allele per row
        :param input_vcf: snpeff annotation results
        :return: tsv file for upload to impala table
        '''
        out_name = "{}.tsv".format(self.get_file_base(input_vcf, 2))
        parse_cmd = 'cat {vcf} | {perl} | java -Xmx16g -jar {snpsift} extractFields \
            - CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].DISTANCE" \
            "ANN[*].HGVS_C" "ANN[*].HGVS_P" > {out}'.format(vcf=input_vcf, perl=self.snpeff_oneperline, \
                                                         snpsift=self.snpsift_jar, out=out_name)
        self.subprocess_cmd(parse_cmd, self.vcf_dir)

    def run_parse(self, vcf_dir):
        '''
        parse snpeff output and create tsv files for
        upload to impala table
        :param vcf_dir: path to directory containing snpeff annotated vcf files
        :return: annotated tsv files for upload to impala
        '''
        for file in os.listdir(vcf_dir):
            if file.endswith('_snpeff.vcf'):
                self.parse_snpeff(file)

# TODO remove chrom prefix

##########  Main Routine  ############
vcf_dir = '/titan/ITMI1/workspaces/users/selasady/impala_scripts/annotation/snpeff'
impala_host = 'glados14'
impala_port = 21050
impala_user_name = 'selasady'

snpeff = snpeff_pipeline(vcf_dir, impala_host, impala_port, impala_user_name)

# snpeff.run_filter(vcf_dir)
# snpeff.run_snpeff(vcf_dir)
snpeff.run_parse(vcf_dir)

