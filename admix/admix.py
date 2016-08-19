#!/usr/bin/env python

"""
Run ADMIXTURE
- create reference marker file
- convert VCF of interest to plink bed format
- run admixture
- output results to impala table
"""

####################
### Script Setup ###
####################

# import modules
import os
import subprocess as sp
import pandas as pd
import sys

class run_admix(object):

    # program paths
    tool_path = '/users/selasady/my_titan_itmi/tools/'
    plink_path = "{}plink1.9/plink".format(tool_path)
    java_path = '/tools/java/jdk1.7/bin/java -jar -Xmx16g'
    gatk_path = "{}GenomeAnalysisTK.jar".format(tool_path)
    admixture_path = "{}admixture_linux-1.3.0/admixture".format(tool_path)
    vcf_verify = "{}snpEff/scripts/vcfBareBones.pl".format(tool_path)
    bcftools = "{}bcftools/bcftools".format(tool_path)
    vcftools = "{}vcftools_0.1.13/bin/vcftools".format(tool_path)

    def __init__(self, vcf_dir, ref_dir, ref_base, panel_file, kval=5, num_cores=15):
        """Returns an admix object with related parameters and files."""
        self.vcf_dir = vcf_dir
        self.ref_dir = ref_dir
        self.ref_base = ref_base
        self.panel_file = panel_file
        self.kval = kval
        self.num_cores = num_cores

    #################################################################################
    ## setup basic methods for working with files and calling subprocess commands  ##
    #################################################################################
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

    @staticmethod
    def create_output_dir(self, input_dir):
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

    ####################################
    ### Create reference marker file ###
    ####################################

    # create panel text file with subject id's, grouped by population
    def create_sorter_and_pop(self):
        '''
        Creates sorted popuplation files for reference data set for super and subpopulations
        :param panel_file: patht to tab-delimited .panel file in format [family_id | subject_id | super_pop | sub_pop ]
        :return sorted_ref.txt a two column tsv ( sample | pop ) to be fed to plink
        for sorting by populations and super_pop.txt and sub_pop.txt
        '''
        out_panel = "{}/sorted_panel.txt".format(self.ref_dir)

        try:
            panel_cols = ['family', 'subject', 'super_pop', 'sub_pop']
            panel = pd.read_csv(self.panel_file, sep='\t', names=panel_cols, index_col=False)
            print panel.head(5)
            # panel.sort_values(by=['super_pop', 'sub_pop'], axis=0, inplace=True)
            # # create super_pop.txt file for super populations
            # super_out = panel[['subject', 'super_pop']]
            # print ("Creating super_pop.txt...")
            # super_out.to_csv("{}/super_pop.txt".format(self.ref_dir), index=False,
            #                  names=['sample', 'super_pop'], sep='\t', header=False)
            # # create sub_pop.txt for subpopulations
            # sub_out = panel[['subject', 'sub_pop']]
            # print ("Creating sub_pop.txt...")
            # sub_out.to_csv("{}/sub_pop.txt".format(self.ref_dir), index=False, names=['sample', 'sub_pop'],
            #                sep='\t', header=False)
            # # create sorted_ref.txt for sorting marker tsv by pop in format family_id | subject_id
            # sorter_out = panel[['family', 'subject']]
            # print("Saving sorted panel file as {}".format(out_panel))
            # sorter_out.to_csv("{}/sorted_panel.txt".format(self.ref_dir), index=False, sep='\t', header=False)
        except Exception as e:
            raise SystemExit("{}".format(e))



#################################################
if __name__ == '__main__':

    admix = run_admix(vcf_dir='/users/selasady/shared/admix/vcf_files/integrated_call_samples_v3.20130502.ALL.panel', \
                      panel_file='/users/selasady/shared/admix/vcf_files/', ref_dir='/users/selasady/shared/admix', \
                      ref_base='all_chroms_1000g_pruned')

    admix.create_sorter_and_pop()



    # # cmd to filter reference ped file to retain MAF between .05 and 0.5 and LD 50 5 0.5
    # self.plink_filter_ped_cmd = "{plink} --file {ref}/{ped} --out {ref}/{ped}_maf_ld --maf 0.05 --max-maf .49 \
    #     --indep-pairwise 50 5 0.5".format(plink=self.plink_path, ref=self.ref_dir, ped=self.ref_base)
    # # cmd to filter reference bed file to retain MAF between .05 and 0.5 and LD 50 5 0.5
    # self.plink_filter_bed_cmd = "{plink} --bfile {ref}/{bed_base} --out {ref}/{bed_base}_maf_ld --maf 0.05 --max-maf .49 \
    #     --indep-pairwise 50 5 0.5".format(plink=self.plink_path, bed_base=self.ref_base, ref=self.ref_dir)
    # # cmd to retain only maf/ld pruned variants and convert ped to bed/bim/fam file using plink
    # self.pruned_file = '{ref}/{ped}_maf_ld.prune.in'.format(ref=self.ref_dir, ped=self.ref_base)
    # # command to convert pruned ped file to bed file
    # self.ped2bed_cmd = '{plink} --file {ref}/{ped} --extract {pruned} --make-bed --indiv-sort file {ref}/sorted_panel.txt \
    #     --out {ref}/{ped}_filtered'.format(plink=run_admix.plink_path, ref=self.ref_dir, ped=self.ref_base,
    #                                        pruned=self.pruned_file)
    # # commend to create filtered bed file from pruned file
    # self.filtered_bed = '{plink} --bfile {ref}/{bed_base} --extract {pruned} --make-bed --indiv-sort file {ref}/sorted_panel.txt \
    #     --out {ref}/{bed_base}_filtered'.format(plink=run_admix.plink_path, ref=self.ref_dir,
    #                                             pruned=self.pruned_file, bed_base=self.ref_base)
    # # create list of snps from marker file to extract from input vcf files
    # self.snp_file = "{}/snplist.txt".format(self.ref_dir)
    # # create output dir for filtered markers
    # self.filtered_out = "{}/filtered/".format(self.vcf_dir)
    # # name filtered bed file output
    # self.marker = "{}/{}.bed".format(self.ref_dir, self.ref_base)