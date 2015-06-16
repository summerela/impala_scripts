#################
# load modules  #
#################
import argparse
import sys
import os
import pandas as pd
from impala.dbapi import connect
from impala.util import as_pandas

#stop irrelevant pandas warnings
pd.options.mode.chained_assignment = None

###################
# parse user args #
###################
try:
    parser = argparse.ArgumentParser()
    parser.add_argument("genes", help="list genes of interest, comma-sep ,no quotes,no spaces", type=str)
    parser.add_argument("db", help="enter database of variants to search", type=str)
    parser.add_argument("kav", help="enter max kaviar frequency to return", type=str)
    parser.add_argument("platform", help="cgi or illumina", type=str)
    args = parser.parse_args()
except:
    e = sys.exc_info()[0]
    print e

##################################################
# read in genes of interest and mark if wildcard #
##################################################
# split gene list by comma
gene_args = args.genes.replace("'", "").split(',')

# mark any wildcards in gene list
wildcards = []
gene_list = []
for gene in gene_args:
    if gene.endswith('%'):
        wildcards.append(gene)
    else:
        gene_list.append(gene)

###############
# Build Query #
###############
print "Building query to look for " + args.platform + " variants in " + str(gene_list) + " and", str(wildcards) + "\n" + \
    "with Kaviar frequency of " + args.kav + "% or less..." + "\n"

# create statements for genes and wildcards
if len(gene_list) > 0 and len(wildcards) < 1:
    gene_statement = "WHERE ens.gene_name IN ('" + "','".join(map(str, gene_list)) + "')"
elif len(wildcards) > 0 and len(gene_list) < 1:
    gene_statement = 'WHERE ens.gene_name LIKE (' + "','".join(map(str, wildcards)) + "')"
elif len(gene_list) > 0 and len(wildcards) > 0:
    gene_statement = "WHERE (ens.gene_name IN ('" + "','".join(
        map(str, gene_list)) + "') OR ens.gene_name LIKE ('" + "'," \
                                                               "'".join(map(str, wildcards)) + "'))"
#build query
if args.platform == 'illumina':
    query = """
            SELECT p.sample_id, p.qual, p.filter, k.id as rsID, (k.alle_freq * 100) as kav_pct, k.alle_cnt as
             kav_count, gene_name, p.chr, p.pos, p.ref, p.alt, p.gt,
                   CASE  WHEN SUBSTRING(p.sample_id, -2) = '01' THEN 'M'
                   WHEN SUBSTRING(p.sample_id, -2) = '02' THEN 'F'
                   WHEN SUBSTRING(p.sample_id, -2) = '03' THEN 'NB' END as member,
                   CONCAT(gene_name, ':', p.chr, ':', CAST(p.pos AS STRING)) as variant_id
             FROM
             (SELECT DISTINCT p.sample_id, p.qual, p.filter, ens.gene_name, p.chr, p.pos, p.ref, p.alt, p.gt
             FROM public_hg19.ensembl_genes ens, {0:s} AS p
             {1:s}
             AND ens.chromosome NOT LIKE 'H%'
             AND ens.chromosome = p.chr
             AND (p.pos >= ens.start AND p.pos <= ens.stop)
             AND p.gt IS NOT NULL  ) AS p
             LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
                   ON p.chr = k.chromosome
                   AND p.pos = k.pos
                   AND p.ref = k.ref
                   AND p.alt = k.alt
             WHERE (k.alle_freq < .{2:s} OR k.alle_freq IS NULL)
            """.format(args.db, gene_statement, args.kav)
elif args.platform == 'cgi':
    query = """
            SELECT p.sample_id, p.allele1varquality, p.totalreadcount, k.id as rsID, (k.alle_freq * 100) as kav_pct,
            k.alle_cnt as kav_count, gene_name, p.chr, p.start, p.stop, p.ref, p.allele1seq, p.allele2seq, p.zygosity,
            (CASE
            WHEN SUBSTRING(p.sample_id, -1) = 'M' THEN 'M'
            WHEN SUBSTRING(p.sample_id, -1) = 'F' THEN 'F'
            WHEN SUBSTRING(p.sample_id, -2) = 'NB' THEN 'NB'
            END) as member, CONCAT(gene_name, ':', p.chr, ':', CAST(p.start AS STRING), ':',
            CAST(p.stop AS STRING)) as variant_id
            FROM
            (SELECT DISTINCT p.sample_id, p.allele1varquality, p.totalreadcount, ens.gene_name, p.chr, \
p.start, p.stop,
             p.ref, p.allele1seq, p.allele2seq, p.zygosity
             FROM public_hg19.ensembl_genes ens, {0:s} AS p
             {1:s}
             AND ens.chromosome NOT LIKE 'H%'
             AND ens.chromosome = p.chr
             AND (p.start >= ens.start AND p.stop <= ens.stop)
             AND p.zygosity IS NOT NULL) AS p
             LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
                ON p.chr = k.chromosome
                AND (k.pos BETWEEN p.start and p.stop)
                AND p.ref = k.ref
                AND p.allele1seq = k.alt
             WHERE (k.alle_freq < .{2:s} OR k.alle_freq IS NULL)
            """.format(args.db, gene_statement, args.kav)
else:
    print "Did you select illumina or cgi as your platform? Please check and try again."

######################
# connect to database #
######################
print "Connecting to impala..." + "\n"

#create database connection
conn = connect(host='glados19', port=21050)
cur = conn.cursor()

#########################
# Run Query on Impala ##
#########################
print "Running query on impala..." + "\n"

#execute query on impala
cur.execute(query)
query_df = as_pandas(cur)

#check that results are returned and print preview
if len(query_df) > 0:
    print "Results found, starting analysis..."
else:
    print "No results found"

############################
# create sets and subsets ##
############################
# create lists to store candidate variants
hom_alt = []
hom_ref = []
comp_het = []
mie = []
missing = []

# group by variant id to locate variants at same chr and pos
by_variant = query_df.groupby('variant_id')

# group variants by gene name
by_gene = query_df.groupby('gene_name')

#######################################
# find hom_alt and hom_ref candidates ##
#######################################
print "Examining homozygous variants from het parents..."

for name, group in by_variant:
    # find positions where both parents are het
    if len(group[(group['member'] == 'M') & (group['gt'] == '0/1')]) > 0 \
        and (len(group[(group['member'] == 'F') & (group['gt'] == '0/1')]) > 0):
        # if the newborn is hom_ref
        if (len(group[(group['member'] == 'NB') & (group['gt'] == '0/0')]) > 0):
            # mark as hom_ref
            group['var_type'] = "hom_ref"
            hom_ref.append(group)
        # if the newborn is hom_alt
        elif (len(group[(group['member'] == 'NB') & (group['gt'] == '1/1')]) > 0):
            group['var_type'] = "hom_alt"
            hom_alt.append(group)
###################
# find comp_hets ##
###################
print "Looking for compound heterozygosity..."

for name, group in by_gene:
    # if there are variants at more than one position and there are nb het variants at diff pos
    if group.pos.nunique() > 1 and (len(group[((group['member'] == 'NB') & (group['gt'] == '0/1'))].pos) > 1):
        # pull out nb het positions
        nb_hets = group[((group['member'] == 'NB') & (group['gt'] == '0/1'))]
        # add parent positions for the gene
        moms = nb_hets.append(group[((group['member'] == "M") & (group['variant_id'].isin(nb_hets['variant_id'])))],
                              ignore_index=True)
        all_vars = moms.append(group[((group['member'] == "F") & (group['variant_id'].isin(nb_hets['variant_id'])))],
                               ignore_index=True)
        # if there are matching variants from both parents, mark as comp-het
        if (len(all_vars[(all_vars['member'] == 'F')]) > 0) and (len(all_vars[(all_vars['member'] == 'M')]) > 0):
            all_vars['var_type'] = "comp_het"
            comp_het.append(all_vars)

#############
# find MIE ##
#############
print "Searching for MIE...."

# het variants
for name, group in by_variant:
    #if newborn is het
    if (len(group[(group['member']=='NB') & (group['gt']=='0/1')]) > 0 ):
        # and both parents at this position are homozygous
        if (len(group[((group['member']=='M') & ((group['gt'] == '0/0')| (group['gt'] == '1/1')))]) > 0) & \
         (len(group[((group['member']=='F') & ((group['gt'] == '0/0')| (group['gt'] == '1/1')))]) > 0):
            if (len(group[(group['member'] == 'F')]) > 0) & (len(group[(group['member'] == 'M')]) > 0):
                group['var_type'] = "mie"
                mie.append(group)
            else:
                group['var_type'] = "missing_parent"
                missing.append(group)
        # or newborn variant is not from either parent
        elif (group[(group['member'] == 'NB')].alt.any() != (group[(group['member'] == 'F')].alt.any() \
                                                                     or group[(group['member'] == 'M')].alt.any())):
            if (len(group[(group['member'] == 'F')]) > 0) & (len(group[(group['member'] == 'M')]) > 0):
                group['var_type'] = "mie"
                mie.append(group)
            else:
                group['var_type'] = "missing_parent"
                missing.append(group)

# hom variants
for name, group in by_variant:
    # if newborn is homozygous
    if (len(group[(group['member']=='NB') & (group['gt']=='1/1')]) > 0 ):
        # and only one parent is het for alt allele
        if ((len(group[( ((group['member']=='F') & (group['gt']=='1/1')) & ((group['member']=='M') & (group['gt']!='1/1')) )]) > 0 ) | (len(group[( ((group['member']=='M') & (group['gt']=='1/1')) & ((group['member']=='F') & (group['gt']!='1/1')) )]) > 0 )):
            if (len(group[(group['member'] == 'F')]) > 0) & (len(group[(group['member'] == 'M')]) > 0):
                group['var_type'] = "mie"
                mie.append(group)
            else:
                group['var_type'] = "missing_parent"
                missing.append(group)
        #or nb variant is not from either parent
        elif (group[(group['member'] == 'NB')].alt.any() != (group[(group['member'] == 'F')].alt.any() \
                                                                     or group[(group['member'] == 'M')].alt.any())):
            if (len(group[(group['member'] == 'F')]) > 0) & (len(group[(group['member'] == 'M')]) > 0):
                group['var_type'] = "mie"
                mie.append(group)
            else:
                group['var_type'] = "missing_parent"
                missing.append(group)

###########################
# merge and save results ##
###########################
print "Saving results to current working directory...."

#create list of groupby result lists and merge
groups = hom_alt + hom_ref + comp_het + mie + missing
results = pd.concat(groups, axis=0, join='outer')

#save to csv
results.to_csv('rare_variant_results.csv', header=True, encoding='utf-8')

#close all the things
cur.close()
sys.exit()

