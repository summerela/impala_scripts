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
import glob
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

    def __init__(self, vcf_dir, ref_dir, ref_base, population='super', kval=5, num_cores=15):
        """Returns an admix object with related parameters and files."""
        self.vcf_dir = vcf_dir
        self.ref_dir = ref_dir
        self.ref_base = ref_base
        self.population = population
        self.kval = kval
        self.num_cores = num_cores
        self.ref_panel = ''.join([f for f in os.listdir(self.ref_dir) if f.endswith('.panel')])
        self.ref_panel_file = "{}/{}".format(self.ref_dir, str(self.ref_panel))
        # cmd to filter reference ped file to retain MAF between .05 and 0.5 and LD 50 5 0.5

        self.plink_filter_bed_cmd = "{plink} --bfile {ref}/{bed_base} --out {ref}/{bed_base}_maf_ld --maf 0.05 --max-maf .49 \
            --indep-pairwise 50 5 0.5".format(plink=self.plink_path, bed_base=self.ref_base, ref=self.ref_dir)
        # cmd to retain only maf/ld pruned variants and convert ped to bed/bim/fam file using plink
        self.pruned_file = '{ref}/{ped}_maf_ld.prune.in'.format(ref=self.ref_dir, ped=self.ref_base)
        self.ped2bed_cmd = '{plink} --file {ref}/{ped} --extract {pruned} --make-bed --indiv-sort file {ref}/sorted_panel.txt \
            --out {ref}/{ped}_filtered'.format(plink=run_admix.plink_path, ref=self.ref_dir, ped=self.ref_base, pruned=self.pruned_file)
        self.filtered_bed = '{plink} --bfile {ref}/{bed_base} --extract {pruned} --make-bed --indiv-sort file {ref}/sorted_panel.txt \
            --out {ref}/{bed_base}_filtered'.format(plink=run_admix.plink_path, ref=self.ref_dir, pruned=self.pruned_file, bed_base=self.ref_base)
        # create list of snps from marker file to extract from input vcf files
        self.snp_file = "{}/snplist.txt".format(self.ref_dir)
        self.filtered_out = "{}/filtered/".format(self.vcf_dir)
        self.marker = "{}/{}.bed".format(self.ref_dir, self.ref_base)

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

    def create_output_dir(self, input_dir):
        if not os.path.exists(self.filtered_out):
            os.makedirs(self.filtered_out)

    ####################################
    ### Create reference marker file ###
    ####################################

    # create text file with subject id's grouped by population for sorting subjects by population
    def create_sorter_and_pop(self, input_dir):
        '''
        ** assumes that reference directory contains a tab-delimited .panel file in format [family_id | subject_id | super_pop | sub_pop ]
        Creates sorted popuplation files for reference data set for super and subpopulations
        :param input_dir: ref_dir containing ref panel file (three columns, tsv, col headers = sample | pop | super_pop )
        :return sorted_ref.txt a two column tsv ( sample | pop ) to be fed to plink
        for sorting by populations and super_pop.txt and sub_pop.txt
        '''
        out_panel = "{}/sorted_panel.txt".format(input_dir)
        panel_check = "{}/*.panel".format(input_dir)
        if glob.glob(panel_check):
            try:
                panel_cols = ['family', 'subject', 'super_pop', 'sub_pop']
                panel = pd.read_csv(self.ref_panel_file, sep='\t', names=panel_cols, index_col=False)
                panel.sort_values(by=['super_pop', 'sub_pop'], axis = 0, inplace=True)
                # create super_pop.txt file for super populations
                super_out = panel[['subject', 'super_pop']]
                print ("Creating super_pop.txt")
                super_out.to_csv("{}/super_pop.txt".format(input_dir), index=False, names=['sample', 'super_pop'], sep='\t', header=False)
                # create sub_pop.txt for subpopulations
                sub_out = panel[['subject', 'sub_pop']]
                print ("Creating sub_pop.txt")
                sub_out.to_csv("{}/sub_pop.txt".format(input_dir), index=False, names=['sample', 'sub_pop'], sep='\t', header=False)
                # create sorted_ref.txt for sorting marker tsv by pop in format family_id | subject_id
                sorter_out = panel[['family', 'subject']]
                print("Saving sorted panel file as {}".format(out_panel))
                sorter_out.to_csv("{}/sorted_panel.txt".format(input_dir), index=False, sep='\t', header=False)
            except Exception as e:
                print e
        else:
            print ("No reference .panel file found.")
            sys.exit(0)

    def ref_snps(self, input_dir):
        # read in reference .bim file
        bim_file = "{ref_dir}/{base}.bim".format(ref_dir=self.ref_dir, base=self.ref_base)
        ref_bim = pd.read_csv(bim_file, sep='\t', header=None, usecols=[0,3], names=['chrom', 'pos'])
        ref_bim.sort_values(by=['chrom', 'pos'], axis=0, inplace=True)
        ref_bim.to_csv("{}/snplist.txt".format(self.ref_dir), header=False, index=False, sep='\t')

    def create_marker_ped(self, input_dir):
        '''
        Runs commands to create marker file from a plink ped file if make_marker == True
        :param input_dir: path to dir containing plink ped/map files
        :return: bed file with MAF and LD filtering to merge with input VCF files for running ADMIXTURE
        '''
        # run the vcf procesing commands
        marker_cmds = [self.plink_filter_ped_cmd,self.ped2bed_cmd]
        for cmd in marker_cmds:
            self.subprocess_cmd(cmd, input_dir)
        self.ref_snps(input_dir)


    def create_marker_bed(self, input_dir):
        '''
        Runs command to create marker file from plink binary bed file if 'make_marker == True'
        :param input_dir: path to dir containing plink bed/bim/fam files
        :return: bed file with MAF and LD Pruning to merge with input vcf files for running ADMIXTURE
        '''
        # marker_cmds = [self.plink_filter_bed_cmd, self.filtered_bed]
        # for cmd in marker_cmds:
        #     self.subprocess_cmd(cmd, input_dir)
        self.ref_snps(input_dir)

    def process_reference(self, input_dir):
        self.create_sorter_and_pop(input_dir)
        bed_check = "{}/*.bed".format(input_dir)
        if glob.glob(bed_check):
            print "\n Found bed file, filtering for LD and MAF. \n"
            self.create_marker_bed(input_dir)
        else:
            print "\n Need to create a bed file. Ensure directory contains ped/map file for reference. \n"
            print("\n Creating marker file... \n")
            self.create_marker_ped(input_dir)

    ###############################
    ### Pre-process VCF file(s) ###
    ###############################

    def make_ref_snps(self, input_dir):
        snp_list = "{}/{}".format(input_dir, ''.join([f for f in os.listdir(input_dir) if f.endswith('.snplist')]))
        snps = pd.read_csv(snp_list, names=['chrom'])
        snp_out = pd.DataFrame(snps['chrom'].apply(lambda x: x.split(':')).values.tolist()).astype(int)
        snp_out.columns = ['chrom', 'pos']
        # add chr prefix to match up 1000g ref with vcf chromosomes
        snp_out['chrom']= 'chr' + snp_out['chrom'].astype(str)
        snp_out.to_csv(self.snp_file, header=None, index=False, sep='\t')

    def subset_vcf(self, input_vcf):
        out_base = self.get_file_base(input_vcf, 2)
        subset_cmd = r'''{vcftools} --gzvcf {vcf_dir}/{input_vcf} --positions {snps} --recode --stdout | {barebones} | gzip > {out_dir}/{out_vcf}_trimmed.vcf.gz'''.format(vcftools=self.vcftools, vcf_dir=self.vcf_dir, input_vcf=input_vcf, snps=self.snp_file, barebones=self.vcf_verify, out_dir=self.filtered_out, out_vcf=out_base)
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
        self.subprocess_cmd(vcf2plink_cmd, self.filtered_out)

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

    def find_beds(self, input_vcf):
        '''
        Locate bed files in a dir and ensure they have matching bim/fam
        :param input_dir: dir where bed files to merge are located
        :return: list of bed files to merge
        '''
        bedList = []
        base_name = self.get_file_base(input_vcf, 2)
        bedList.append(str(self.vcf_dir + '/' + base_name + '.bed'))
        bedList.append(str(self.vcf_dir + '/' + base_name + '.bim'))
        bedList.append(str(self.vcf_dir + '/' + base_name + '.fam'))
        bedList.sort()
        return bedList

    def plink2_merge(self, input_vcf):
        '''
        Generate command to merge bed files using plink then run with subprocess
        :param input_vcf:
        :return: merged bed/bim/fam file of 1000g marker plus input bed files
        '''
        print ("\n Merging bed file(s) with marker file... \n")
        base_name = self.get_file_base(input_vcf, 2)
        beds = self.find_beds(input_vcf)
        bed_list = ' '.join(beds)
        merge_cmd = "{} --bfile {}/{} --bmerge {} --memory 40000 --out {}/{}_merged".format(self.plink_path, self.ref_dir, self.ref_base, bed_list, self.filtered_out, base_name)
        self.subprocess_cmd(merge_cmd, self.filtered_out)


    def run_merge(self, input_dir):
        for file in os.listdir(input_dir):
            if file.endswith('_trimmed.vcf.gz'):
                bed_list = self.find_beds(file)
                self.plink2_merge(file)

    #######################
    ###  run admixture  ###
    #######################

    def make_pop_file(self, input_bed, kval):
        '''
        Create .pop file with 5 known major populations from 1000g and '-' for unknown from vcf files
        Assumes ref_dir contains tsv files in format subject_id | pop named super_pop.txt and sub_pop.txt
        :param input_vcf: input vcf file to grab base file name
        :param kval: specify 5 for super population and 26 for subpopulations
        :return .pop file for merged bed/bim/fam set as input for supervised admixture
        '''
        # read in _merged.fam file
        base_name =self.get_file_base(input_bed, 1)
        in_fam = "{}/{}.fam".format(self.filtered_out, base_name)
        merged_df = pd.read_csv(in_fam, sep='\t', usecols=[1], names=['sample'])
        if self.population == 'super':
            # merge .fam file with pop file to get population for each known subject
            super_in = "{}/super_pop.txt".format(self.ref_dir)
            super_map = pd.read_csv(super_in, sep='\t', header=None, names=['sample', 'super_pop'])
            super_out = pd.merge(merged_df, super_map, how='left', on = 'sample')
            super_out.drop('sample', inplace=True, axis=1)
            super_file = "{}/{}.pop".format(self.filtered_out, base_name)
            print ("Saving pop file as {}".format(super_file))
            super_out.to_csv(super_file, sep='\t', na_rep='-', header=False, index=False)
        elif self.population == 'sub':
            # create sub_pop file
            sub_in = "{}/sub_pop.txt".format(self.ref_dir)
            sub_map = pd.read_csv(sub_in, sep='\t', header=None, names=['sample', 'sub_pop'])
            sub_out = pd.merge(merged_df, sub_map, how='left', on= 'sample')
            sub_out.drop('sample', inplace=True, axis=1)
            sub_file = "{}/{}.pop".format(self.filtered_out, base_name)
            print ("Saving pop file as {}".format(sub_file))
            sub_out.to_csv(sub_file, sep='\t', na_rep='-', header=False, index=False)
        else:
            print ("Please select either super or sub for your admixture population.")

    def make_pop(self, input_dir, kval):
        for file in os.listdir(input_dir):
            if file.endswith('_merged.bed'):
                self.make_pop_file(file, kval)

    def admix_cmd(self, num_cores, merged_bed, kval):
        '''
        run admixture on input merged bed files
        :param num_cores: number of cores to use for running admixture
        :param merged_bed: bed file of merged marker plus unknown
        :param kval: choose either 5 or 26 to run admix with super or sub populations
        :return: P and Q files with admixture results
        '''
        admix_cmd = "nohup {} -j{} {} --supervised {}".format(self.admixture_path, int(num_cores), merged_bed, kval)
        self.subprocess_cmd(admix_cmd, self.filtered_out)

    def run_admix_program(self, input_dir, num_cores, kval):
        for file in os.listdir(input_dir):
            if file.endswith('_merged.bed'):
                self.admix_cmd(num_cores, file, kval)

    ######################
    ## Process Results ###
    ######################

    def get_results(self, q_file):
        # read admixture results in
        q_file_in = pd.read_csv("{}/{}".format(admix.filtered_out, q_file), header=None, sep=' ')
        pop_in = "{}/{}.pop".format(admix.filtered_out, admix.get_file_base(q_file, 2))
        # read in pop file
        pop_file = pd.read_csv(pop_in, header=None, sep=' ', names=['pop'])
        # read in fam file
        fam_file_in = "{}/{}.fam".format(admix.filtered_out, admix.get_file_base(q_file, 2))
        fam_file = pd.read_csv(fam_file_in, header=None, sep='\t', usecols=[1], names=['sample'])
        frames = [fam_file, pop_file, q_file_in]
        out_frame = pd.concat(frames, ignore_index=True, keys=None, axis=1)
        pop_list = list(out_frame[1][(out_frame[1] != '-')].unique())
        head_list = ['sample', 'known_pop']
        col_list = head_list + pop_list
        out_frame.columns = col_list
        results_df = out_frame[(out_frame['known_pop'] == '-')]
        results_df.drop(['known_pop'],inplace=True,axis=1,errors='ignore')
        return results_df

    def create_results_file(self, input_dir):
        for file in os.listdir(input_dir):
            if file.endswith('.Q'):
                result = self.get_results(file)
                result.to_csv("{}/results.txt".format(admix.filtered_out), header=True, sep='\t', index=False, mode='a')

########################################################

if __name__ == '__main__':

    ## edit the file paths below for each run ##

    # path to input vcf files
    vcf_file_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/test4/'
    # path to reference ped/map files
    # ref_file_dir = '/users/selasady/my_titan_itmi/impala_scripts/annotation/admix/1000g_vcf/'
    ref_file_dir = '/titan/ITMI1/workspaces/users/selasady/impala_scripts/annotation/admix/blog_ref/'
    # ped/map file basename
    # ref_base_name = 'ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05'
    ref_base_name = 'ref'
    # choose to run admix for either the super or sub populations found in your reference set
    population = 'super'
    # specify the number of distinct populations in your chosen reference set
    pop_kval = 17
    # number of cores to run admixture
    admix_cores = 20
    # mark true if reference marker needs to be created, else False
    create_marker = 'True'

    ##########  Main Routine  ############

    admix = run_admix(vcf_file_dir, ref_file_dir, ref_base_name, population, pop_kval, admix_cores)

    def run_admixture(cores, kval):
        admix.filter_vcf(admix.vcf_dir)
        admix.convert_vcf(admix.filtered_out)
        admix.run_merge(admix.filtered_out)
        admix.make_pop(admix.filtered_out, admix.population)
        admix.run_admix_program(admix.filtered_out, admix_cores, pop_kval)
        admix.create_results_file(admix.filtered_out)

    ### create marker ###
    if create_marker == 'True':
        admix.create_output_dir(admix.vcf_dir)
        admix.process_reference(admix.ref_dir)
        run_admixture(admix_cores, population)
    else:
        admix.create_output_dir(admix.vcf_dir)
        run_admixture(admix_cores, population)



