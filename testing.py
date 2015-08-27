__author__ = 'selasady'

import pandas as pd

# import txt file from ucsc
kg = pd.read_table('/Users/selasady/Downloads/knownGene.txt')

# name columns
kg.columns =['ucsc_id', 'chrom', 'strand', 'tx_start', 'tx_stop', 'cds_start','cds_stop', 'exon_count','exon_starts', 'exon_stops', 'protein_id', 'align_id']

# remove dangling commas from ends of lists
#kg['exon_starts'] = kg['exon_starts'].str.rstrip(',')
kg['exon_stops'] = kg['exon_stops'].str.rstrip(',')

# remove 'chr' prefix from chromosome field
kg['chrom']= kg['chrom'].str.replace('chr', '')

# add 1 to start coord to convert from half-open zero to one-based
kg['tx_start'] = kg['tx_start'] + 1

# add 1 to each exon start site
new_exon_starts = []
for x in kg['exon_starts']:
    value_list = []
    for y in x.rstrip(',').split(','):
        value_list.append(str(int(y) + 1))
    new_exon_starts.append(','.join(value_list))
kg['exon_starts'] = new_exon_starts

# # save as tsv file
kg.to_csv('ucsc_knowngene.tsv', index=False, sep='\t')
