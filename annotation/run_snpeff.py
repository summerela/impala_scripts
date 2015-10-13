'''
run_snpeff.py
purpose: add functional annotation to variants in vcf format using SnpEff
input: vcf file
    user arguments:
                   snpeff_path = path to snpeff.jar file
                   snpeff_ref = reference for annotation(ex.GRCh37.75)
                   vcf_outname = desired name for output vcf file
output: tab delimited file with added functional annotations
'''

# create command to run
snpeff_cmd = 'java -Xmx4g -jar {} -v {} {} > nbs_snpeff.vcf'.format(snpeff_path, snpeff_ref, vcf_outname)
# run snpeff on vcf file
snpeff_process = subp.Popen(snpeff_cmd, shell=True)
snpeff_process.wait()