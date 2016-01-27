###################
## set variables ##
###################

# specify file name prefix
in_name = 'snpeff_vars'

# set working directory to location of vcf files
import os
os.chdir('/users/selasady/my_titan_itmi/impala_scripts/testing/')

# set file paths to exectuables
plink_path = '/users/selasady/my_titan_itmi/tools/plink/plink --noweb'
tabix_path = '/users/selasady/my_titan_itmi/tools/tabix-0.2.6/tabix'

##########################
## convert vcf to plink ##
##########################

# convert vcf to plink format with vcftools
for file in os.listdir(os.getcwd()):
    if file.endswith(in_name + '.vcf'):
        print "Converting {} to plink format... \n".format(file)
        plink_out = str('.'.join(file.split('.')[:-1]) if '.' in file else file) + '.plink'
        vcf2plink_cmd = 'cat {} | {} > {}'.format(file, vcf_verify,vcf_checked_out)
        # create the file and run snpeff
        with open(vcf_checked_out, "w") as out_file:
            try:
                ps = subprocess.Popen(snp_verify_cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                print ps.communicate()[0]
            except subprocess.CalledProcessError as e:
                 print e.output



for file in *.vcf.gz; \
         do /Users/summerrae/tools/plink --vcf $file -recode 12 --out ${f%%.*} --memory 30000 --threads 8 --no-fid --no-parents --no-sex --no-pheno --double-id --geno 0.1; done