__author__ = 'Summer Elasady for ISB'

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
    parser = argparse.ArgumentParser(description='Find pathogenic variants using clinvar.', usage='%(prog)s [options]')
    parser.add_argument("--chr", help="chromosome as 7,8,9 or 'all'", type=str, nargs='?', default='all')
    parser.add_argument("--member", help="trio members as M,F,NB; or 'all'", type=str, nargs='?', default='all')
    parser.add_argument("--sample_id", help="sample id's as 101-454-M,101-345-F or 'all'", type=str, nargs='?',
                        default='all')
    parser.add_argument("--kav", help="max kaviar frequency as 10", type=str, nargs='?', default='10')
    parser.add_argument("--platform", help="cgi or illumina", type=str, nargs='?', default='illumina')
    args = parser.parse_args()
except:
    e = sys.exc_info()[0]
    print e

####################
# format user args #
####################

def process_args(arg, val):
    if val != 'all' and (',' in val):
        print "WHERE clin.{0} IN ('" + "'".join(map(str, val)) + "')".format(arg)
    else:
        print "WHERE clin.{0} = {1}".format(arg,val)

for key, value in vars(args).items():
    process_args(key, value)




# gene_args = args.genes.replace("'", "").split(',')
#
# # mark any wildcards in gene list
# wildcards = []
# gene_list = []
# for gene in gene_args:
#     if gene.endswith('%'):
#         wildcards.append(gene)
#     else:
#         gene_list.append(gene)

# ###############
# # Build Query #
# ###############
# print "Building query to look for " + args.platform + " variants in " + str(gene_list) + " and", str(wildcards) + "\n" + \
#     "with Kaviar frequency of " + args.kav + "% or less..." + "\n"
#
# # create statements for genes and wildcards
# if len(gene_list) > 0 and len(wildcards) < 1:
#     gene_statement = "WHERE ens.gene_name IN ('" + "','".join(map(str, gene_list)) + "')"
# elif len(wildcards) > 0 and len(gene_list) < 1:
#     gene_statement = 'WHERE ens.gene_name LIKE (' + "','".join(map(str, wildcards)) + "')"
# elif len(gene_list) > 0 and len(wildcards) > 0:
#     gene_statement = "WHERE (ens.gene_name IN ('" + "','".join(
#         map(str, gene_list)) + "') OR ens.gene_name LIKE ('" + "'," \
#                                                                "'".join(map(str, wildcards)) + "'))"
# #build query
# if args.platform == 'illumina':
#     query = """
#             SELECT p.sample_id, p.qual, p.filter, k.id as rsID, (k.alle_freq * 100) as kav_pct, k.alle_cnt as
#              kav_count, gene_name, p.chr, p.pos, p.ref, p.alt, p.gt,
#                    CASE  WHEN SUBSTRING(p.sample_id, -2) = '01' THEN 'M'
#                    WHEN SUBSTRING(p.sample_id, -2) = '02' THEN 'F'
#                    WHEN SUBSTRING(p.sample_id, -2) = '03' THEN 'NB' END as member,
#                    CONCAT(gene_name, ':', p.chr, ':', CAST(p.pos AS STRING)) as variant_id
#              FROM
#              (SELECT DISTINCT p.sample_id, p.qual, p.filter, ens.gene_name, p.chr, p.pos, p.ref, p.alt, p.gt
#              FROM public_hg19.ensembl_genes ens, {0:s} AS p
#              {1:s}
#              AND ens.chromosome NOT LIKE 'H%'
#              AND ens.chromosome = p.chr
#              AND (p.pos >= ens.start AND p.pos <= ens.stop)
#              AND p.gt IS NOT NULL  ) AS p
#              LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
#                    ON p.chr = k.chromosome
#                    AND p.pos = k.pos
#                    AND p.ref = k.ref
#                    AND p.alt = k.alt
#              WHERE (k.alle_freq < .{2:s} OR k.alle_freq IS NULL)
#             """.format(args.db, gene_statement, args.kav)
# elif args.platform == 'cgi':
#     query = """
#             SELECT p.sample_id, p.allele1varquality, p.totalreadcount, k.id as rsID, (k.alle_freq * 100) as kav_pct,
#             k.alle_cnt as kav_count, gene_name, p.chr, p.start, p.stop, p.ref, p.allele1seq, p.allele2seq, p.zygosity,
#             (CASE
#             WHEN SUBSTRING(p.sample_id, -1) = 'M' THEN 'M'
#             WHEN SUBSTRING(p.sample_id, -1) = 'F' THEN 'F'
#             WHEN SUBSTRING(p.sample_id, -2) = 'NB' THEN 'NB'
#             END) as member, CONCAT(gene_name, ':', p.chr, ':', CAST(p.start AS STRING), ':',
#             CAST(p.stop AS STRING)) as variant_id
#             FROM
#             (SELECT DISTINCT p.sample_id, p.allele1varquality, p.totalreadcount, ens.gene_name, p.chr, \
# p.start, p.stop,
#              p.ref, p.allele1seq, p.allele2seq, p.zygosity
#              FROM public_hg19.ensembl_genes ens, {0:s} AS p
#              {1:s}
#              AND ens.chromosome NOT LIKE 'H%'
#              AND ens.chromosome = p.chr
#              AND (p.start >= ens.start AND p.stop <= ens.stop)
#              AND p.zygosity IS NOT NULL) AS p
#              LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
#                 ON p.chr = k.chromosome
#                 AND (k.pos BETWEEN p.start and p.stop)
#                 AND p.ref = k.ref
#                 AND p.allele1seq = k.alt
#              WHERE (k.alle_freq < .{2:s} OR k.alle_freq IS NULL)
#             """.format(args.db, gene_statement, args.kav)
# else:
#     print "Did you select illumina or cgi as your platform? Please check and try again."
#
# ######################
# # connect to database #
# ######################
# print "Connecting to impala..." + "\n"
#
# #create database connection
# conn = connect(host='glados19', port=21050)
# cur = conn.cursor()
#
# #########################
# # Run Query on Impala ##
# #########################
# print "Running query on impala..." + "\n"
#
# #execute query on impala
# cur.execute(query)
# query_df = as_pandas(cur)
#
# #check that results are returned and print preview
# if len(query_df) > 0:
#     print "Results found, starting analysis..."
# else:
#     print "No results found"
#
