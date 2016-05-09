#!/usr/bin/env python

import pandas as pd
import csv

####################
## 23andme to VCF ##
####################

# read in file
in_file = "D:\Documents\julian_genome.txt"

# create header
with open('D:\Documents\julian.vcf', 'w') as csvfile:
    csvfile.write("##fileformat=VCFv4.0")
    csvfile.write('\n')
    csvfile.write("##Patient=Julian")
    csvfile.write('\n')

in_df = pd.read_csv(in_file, sep='\t', comment='#', names=['ID', '#CHROM', 'POS', 'GENO'])
out_df = in_df[['#CHROM', 'POS', 'ID']]
out_df['REF'] = in_df['GENO'].str[0]
out_df['ALT'] = in_df['GENO'].str[1:]
out_df['QUAL'] = '100'
out_df['FILTER'] = 'PASS'
out_df['INFO'] = 'AC=2;AF=0.00122;AN=1644;DS;set=Intersection'

# out_df.to_csv('D:\Documents\julian.vcf', sep='\t', mode='a', index=False)

## filter results ##

results_in = pd.read_csv('D:\Documents\julian_genome_results.csv', sep=',')




