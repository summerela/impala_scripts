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
import csv
import datetime


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

    def __init__(self, vcf_dir, ref_dir, ped_base, kval=5, num_cores=15):
        """Returns an admix object with related parameters and files."""
        self.vcf_dir = vcf_dir
        self.ref_dir = ref_dir
        self.ped_base = ped_base
        self.kval = kval
        self.num_cores = num_cores
        self.ref_panel = ''.join([f for f in os.listdir(self.ref_dir) if f.endswith('.panel')])
        self.ref_panel_file = "{}/{}".format(self.ref_dir, str(self.ref_panel))
        # cmd to filter reference ped file to retain MAF between .05 and 0.5 and LD 50 5 0.5
        self.plink_filter_cmd = "{plink} --file {ref}/{ped} --out {ref}/{ped}_maf_ld --maf 0.05 --max-maf .49 \
            --indep-pairwise 50 5 0.5".format(plink=run_admix.plink_path, ref=ref_dir, ped=self.ped_base)
        # cmd to retain only maf/ld pruned variants and convert ped to bed/bim/fam file using plink
        self.pruned_file = '{ref}/{ped}_maf_ld.prune.in'.format(ref=self.ref_dir, ped=self.ped_base)
        self.ped2bed_cmd = '{plink} --file {ref}/{ped} --extract {pruned} --make-bed --indiv-sort file {ref}/sorted_panel.txt \
            --out {ref}/{ped}'.format(plink=run_admix.plink_path, ref=self.ref_dir, ped=self.ped_base, pruned=self.pruned_file)
        # create list of snps from marker file to extract from input vcf files
        self.ref_snps = "{plink} --bfile {ped} --write-snplist --out {ped}_snps".format(plink=run_admix.plink_path, ped=self.ped_base)
        self.snp_file = "{}/snplist.txt".format(self.ref_dir)
        self.filtered_out = "{}/filtered/".format(self.vcf_dir)
        self.marker = "{}/{}.bed".format(self.ref_dir, self.ped_base)

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
        ps = sp.Popen(command, shell=True,stdout=sp.PIPE,stderr=sp.PIPE, cwd=input_dir)
        try:
           print ps.communicate()
        except sp.CalledProcessError as e:
             print e

    @staticmethod
    def get_file_base(in_file, num_to_remove):
        basename = str('.'.join(in_file.split('.')[:-int(num_to_remove)]) if '.' in in_file else in_file)
        return basename

    ####################################
    ### Create reference marker file ###
    ####################################

    # create text file with subject id's grouped by population for sorting subjects by population
    def create_sorter_and_pop(self, input_dir):
        '''
        :param input_dir: ref_dir containing ref panel file (three columns, tsv, col headers = sample | pop | super_pop )
        :return sorted_ref.txt a two column tsv ( sample | pop ) to be fed to plink for sorting by pop
        '''
        out_panel = "{}/sorted_panel.txt".format(input_dir)
        print("Saving sorted panel file as {}".format(out_panel))
        try:
            panel = pd.read_csv(self.ref_panel_file, sep='\t', usecols=['sample', 'pop', 'super_pop', 'gender'])
            panel.sort_values(by=['super_pop', 'pop'], axis = 0, inplace=True)
            # create super_pop.txt file for super populations
            super_out = panel[['sample', 'pop']]
            super_out.to_csv("{}/super_pop.txt".format(input_dir), index=False, names=['sample', 'super_pop'])
            return super_out
            # create sub_pop.txt for subpopulations
            sub_out = panel[['sample', 'super_pop']]
            sub_out.to_csv("{}/sub_pop.txt".format(input_dir), index=False, names=['sample', 'pop'])
            # create sorted_ref.txt for sorting marker tsv by pop in format family_id | subject_id
            sorter_out = panel[['sample', 'sample']]
            sorter_out.to_csv("{}/sorted_panel.txt".format(input_dir), index=False, sep='\t', header=False)
        except Exception as e:
            print e

    def create_marker_file(self, input_dir):
        '''
        Runs commands to create marker file if make_marker == True
        :param plink: path to plink1.9
        :param ref_dir: path to dir containing plink ped/map files
        :return: marker file to merge with input VCF files for running ADMIXTURE
        '''
        # run the vcf procesing commands
        marker_cmds = [self.plink_filter_cmd,self.ped2bed_cmd, self.ref_snps]
        for cmd in marker_cmds:
            self.subprocess_cmd(cmd, input_dir)

    def process_reference(self, input_dir):
        print("\n Creating sorted panel file and pop files for super and sub populations... \n")
        self.create_sorter_and_pop(input_dir)
        print("\n Creating marker file... \n")
        self.create_marker_file(input_dir)

    ###############################
    ### Pre-process VCF file(s) ###
    ###############################

    def make_ref_snps(self, input_dir):
        snp_list = "{}/{}".format(input_dir, ''.join([f for f in os.listdir(input_dir) if f.endswith('.snplist')]))
        snps = pd.read_csv(snp_list, names=['chrom'])
        snp_out = pd.DataFrame(snps['chrom'].apply(lambda x: x.split(':')).values.tolist()).astype(int)
        snp_out.columns = ['chrom', 'pos']
        snp_out['chrom']= 'chr' + snp_out['chrom'].astype(str)
        snp_out.to_csv(self.snp_file, header=None, index=False, sep='\t')

    def subset_vcf(self, input_vcf):
        out_base = self.get_file_base(input_vcf, 2)
        subset_cmd = r'''{vcftools} --gzvcf {vcf_dir}/{input_vcf} --positions {snps} --recode --stdout | {barebones} | gzip > {out_vcf}_trimmed.vcf.gz'''.format(vcftools=self.vcftools, vcf_dir=self.vcf_dir, input_vcf=input_vcf, snps=self.snp_file, barebones=self.vcf_verify, out_vcf=out_base)
        self.subprocess_cmd(subset_cmd, self.vcf_dir)

    def vcf_to_bed(self, input_vcf):
        '''
        convert stripped vcf to bed/bim/fam with plink 1.9
        :param input_vcf: vcf file to process
        :return: stripped, verified vcf files converted to plink binary (bed) format
        '''
        bed_out = self.get_file_base(input_vcf, 2)
        vcf2plink_cmd = "nohup {plink} --vcf {filtered_dir}/{file} --double-id --biallelic-only strict --geno 0.1 --allow-no-sex --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out {filtered_dir}/{bed}".format(plink=self.plink_path, \
                                        ref=self.ref_dir, snp=self.snp_file, filtered_dir=self.filtered_out, file=input_vcf, bed=bed_out)
        print ("Converting {} from vcf to bed/bim/fam").format(input_vcf)
        # self.subprocess_cmd(vcf2plink_cmd, self.filtered_out)
        print vcf2plink_cmd


    def filter_vcf(self, vcf_dir):
        for file in os.listdir(vcf_dir):
            if file.endswith('.vcf.gz'):
                self.subset_vcf(file)

    def convert_vcf(self, filtered_dir):
        for file in os.listdir(filtered_dir):
            if file.endswith('_trimmed.vcf.gz'):
                self.vcf_to_bed(file)

    #####################################
    ### Merge BED file(s) with Marker ###
    #####################################

    def create_merge_list(self, input_dir):
        '''
        Locate bed files in a dir and ensure they have matching bim/fam
        :param input_dir: dir where bed files to merge are located
        :return: list of bed files to merge
        '''
        bedList = []
        for file in os.listdir(input_dir):
            if file.endswith('_trimmed.vcf.gz'):
                base_name = "{}.bed".format(self.get_file_base(file, 2))
                bedList.append(base_name)
                bedList.sort()
        plinkList = []
        for bed in bedList:
            (stem, ext) = os.path.splitext(bed)
            plinkList.append(stem)
            for suffix in (".bim", ".fam"):
                myFile = stem+suffix
                if not os.path.exists(input_dir + '/' + myFile):
                    print ("Missing Plink data file "+myFile)
        return plinkList

    def create_merge_text(self, input_vcf, merge_list):
        '''
        create a three column text file of bed/bim/fam files to merge based
        on list from findBedStems
        :param input_vcf: vcf file to merge with marker file
        :param merge_list: list generated with findBedStems()
        :return: tsv file of each bed bim fam file to  merge for input to plink merge
        '''
        print ("Creating merge file for {}".format(input_vcf))
        file_list = []
        for item in merge_list:
            bed_out = str(item + '.bed')
            bim_out = str(item + '.bim')
            fam_out = str(item + '.fam')
            ref_bed = "{}/{}.bed".format(self.ref_dir, self.ped_base)
            ref_bim = "{}/{}.bim".format(self.ref_dir, self.ped_base)
            ref_fam = "{}/{}.fam".format(self.ref_dir, self.ped_base)
            t = [bed_out, bim_out, fam_out]
            r = [ref_bed, ref_bim, ref_fam]
            file_list.append(t)
            file_list.append(r)
        vcf_base = self.get_file_base(input_vcf, 2)
        out_file = "{}/{}_merge.txt".format(self.filtered_out, vcf_base)
        with open(out_file,'wb') as outf:
            writer = csv.writer(outf, delimiter='\t')
            writer.writerows(file_list)

     def plink2_merge(self, input_vcf):
        '''
        Generate command to merge bed files using plink then run with subprocess
        :param input_vcf:
        :return: merged bed/bim/fam file of 1000g marker plus input bed files
        '''
        print ("Merging bed file(s) with marker file... \n")
        vcf_base = self.get_file_base(input_vcf, 2)
        merge_file = "{}/{}_merge.txt".format(self.filtered_out, vcf_base)
        merge_cmd = "{} --merge-list {} --make-bed --memory 40000 --merge-equal-pos --out {}/{}_merged".format(self.plink_path, merge_file, self.filtered_dir, vcf_base)
        self.subprocess_cmd(merge_cmd, self.filtered_dir)


    def setup_merge(self, input_dir):
        bed_list = self.create_merge_list(self.filtered_out)
        for file in os.listdir(input_dir):
            if file.endswith('_trimmed.vcf.gz'):
                self.create_merge_text(file, bed_list)

    def merge_beds(self, filtered_dir):
        for file in os.listdir(filtered_dir):






########################################################

if __name__ == '__main__':
    ### File paths ###

    # path to input vcf files
    vcf_file_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test3/'
    # path to reference ped/map files
    ref_file_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'
    # ped/map file basename
    ped_base_name = 'ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
    # choose either 5 to run admix with super pops, 26 for sub-populations, or 'both'
    pop_kval = 5
    # number of cores to run admixture
    admix_cores = 20

##########  testing   ############

admix = run_admix(vcf_file_dir, ref_file_dir, ped_base_name, pop_kval, admix_cores)

# run stuff
# admix.process_reference(admix.ref_dir)
# admix.make_ref_snps(admix.ref_dir)
# admix.filter_vcf(admix.vcf_dir)
# admix.convert_vcf(admix.filtered_out)
# admix.merge_beds(admix.filtered_out)
admix.plink2_merge(admix.filtered_out)











    # #######################
    # ###  run admixture  ###
    # #######################
    #
    # def make_pop(ref_dir, out_dir, merged_fam, kval):
    #     '''
    #     Create .pop file with 5 known major populations from 1000g and '-' for unknown from vcf files
    #     Assumes ref_dir contains tsv files in format subject_id | pop named super_pop.txt and sub_pop.txt
    #     :param out_dir path to directory containing merged fam files
    #     :param ref_dir: path to location of reference files used to create marker set
    #     :param merged_fam: path to merged fam file to create pop file for
    #     :param kval: specify 5 for super population and 26 for subpopulations
    #     :return .pop file for merged bed/bim/fam set as input for supervised admixture
    #     '''
    #     # read in _merged.fam file
    #     in_fam = "{}/{}".format(out_dir, merged_fam)
    #     merged_df = pd.read_csv(in_fam, sep=' ', usecols=[1], names=['sample'])
    #     fam_basename = str('.'.join(merged_fam.split('.')[:-1]) if '.' in merged_fam else merged_fam)
    #     if kval == 5:
    #         # pull population from the 1000g map file
    #         in_file = "{}/super_pop.txt".format(ref_dir)
    #         ref_map = pd.read_csv(in_file, sep='\t', header=0)
    #         out_df = pd.merge(merged_df, ref_map, how='left', on = 'sample')
    #         out_df = out_df['super_pop']
    #         out_file = "{}/{}.pop".format(out_dir, fam_basename)
    #         print ("Saving pop file as {}".format(out_file))
    #         out_df.to_csv(out_file, sep='\t', na_rep='-', header=False, index=False)
    #     elif kval == 26:
    #         # pull population from the 1000g map file
    #         in_file = "{}/sub_pop.txt".format(ref_dir)
    #         ref_map = pd.read_csv(in_file, sep='\t', header=0)
    #         out_df = pd.merge(merged_df, ref_map, how='left', left_on= 'Individual ID', right_on= 'sample')
    #         out_df = out_df['pop']
    #         out_file = "{}/{}.pop".format(out_dir, fam_basename)
    #         print ("Saving pop file as {}".format(out_file))
    #         out_df.to_csv(out_file, sep='\t', na_rep='-', header=False, index=False)
    #
    # def admix_cmd(admix_path, num_cores, merged_bed, out_dir, kval):
    #     '''
    #     run admixture on input merged bed files
    #     :param admix_path: path to admixture program
    #     :param num_cores: number of cores to use for running admixture
    #     :param merged_bed: bed file of merged marker plus unknown
    #     :param out_dir: path to bed file directory
    #     :param kval: choose either 5 or 26 to run admix with super or sub populations
    #     :return: P and Q files with admixture results
    #     '''
    #     admix_cmd = "nohup {} -j{} {} --supervised {}".format(admixture_path, int(num_cores), merged_bed, kval)
    #     subprocess_cmd(admix_cmd, out_dir)
    #
    # def run_pop(ref_dir, out_dir, kval):
    #     for file in os.listdir(out_dir):
    #         if file.endswith('_merged.fam'):
    #             make_pop(ref_dir, out_dir, file, kval)
    #
    # def run_admix(admix_path, num_cores, out_dir, kval):
    #     for file in os.listdir(out_dir):
    #         if file.endswith('_merged.bed'):
    #             admix_cmd(admix_path, num_cores, file, out_dir, kval)

    #######################
    ### Process Results ###
    #######################

    # # mark true if reference marker needs to be created, else False
    # create_marker = 'False'
    #
    # ### Run ADMIXTURE ###
    #
    # make_out_dir(vcf_file_dir)
    # #filtered_out = "{}/filtered/".format(vcf_file_dir)
    #
    # if create_marker == 'True':
    #     print ("Making marker... ")
    #     process_reference(ref_file_dir, plink_path)
    #
    # process_vcf(vcf_file_dir, vcf_verify, ref_file_dir)

    # for file in os.listdir(filtered_out):
    #     if file.endswith(".Q"):
    #         # read admixture results in
    #         q_file = pd.read_csv("{}/{}".format(filtered_out, file), header=None, sep=' ')
    #         base_start = str('.'.join(file.split('.')[:-2]) if '.' in file else file)
    #         pop_file = "{}/{}.pop".format(filtered_out, base_start)
    #         # read in pop file
    #         pop_file = pd.read_csv(pop_file, header=None, sep=' ')
    #         # read in fam file
    #         fam_file = "{}/{}.fam".format(filtered_out, base_start)
    #         fam_file = pd.read_csv(fam_file, header=None, sep= ' ', usecols=[1])
    #         frames = [fam_file, pop_file, q_file]
    #         out_frame = pd.concat(frames, ignore_index=True, keys=None, axis=1)
    #         out_frame.columns= ['subject', 'known_pop', 'EUR', 'EAS', 'AMR', 'SAS', 'AFR']
    #         print out_frame[(out_frame['known_pop'] == '-')]
