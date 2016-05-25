#!/usr/bin/env python



vcf_dir = '/titan/ITMI1/workspaces/users/selasady/impala_scripts/annotation/snpeff'
impala_host = 'glados14'
impala_port = 21050
impala_user_name = 'selasady'
hdfs_path = '/user/selasady/'

snpeff = snpeff_pipeline(vcf_dir, impala_host, impala_port, impala_user_name, hdfs_path)

#######################
# run snpeff routines #
#######################
snpeff.run_snpeff_routine()
snpeff.cur.close()