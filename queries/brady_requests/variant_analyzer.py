# _author__ = 'summerrae'
# ######################
# # process user args  #
# ######################
# import argparse
# import sys
#
# try:
#     parser = argparse.ArgumentParser()
#     parser.add_argument("genes", help="list genes of interest, comma-sep ,no quotes,no spaces", type=str)
#     parser.add_argument("db", help="enter database of variants to search", type=str)
#     parser.add_argument("kav", help="enter max kaviar frequency to return", type=str)
#     parser.add_argument("platform", help="cgi or illumina", type=str)
#     args = parser.parse_args()
# except:
#     e = sys.exc_info()[0]
#     print e
#
# ##################################################
# # read in genes of interest and mark if wildcard #
# ##################################################
# # split gene list by comma
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
#
# ###############
# # Build Query #
# ###############
# print "Building query to look for " + args.platform + " variants in " + str(gene_list) + "\n" + "with Kaviar frequency of " + args.kav + "% or less..." + "\n"
#
# # create statements for genes and wildcards
# if len(gene_list) > 0 and len(wildcards) < 1:
#     gene_statement = "WHERE ens.gene_name IN ('" + "','".join(map(str, gene_list)) + "')"
# if len(wildcards) > 0 and len(gene_list) < 1:
#     gene_statement = 'WHERE ens.gene_name LIKE (' + "','".join(map(str, wildcards)) + "')"
# if len(gene_list) > 0 and len(wildcards) > 0:
#     gene_statement = "WHERE (ens.gene_name IN ('" + "','".join(
#         map(str, gene_list)) + "') OR ens.gene_name LIKE ('" + "'," \
#                                                                "'".join(map(str, wildcards)) + "'))"
# #build illumina query
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
#              WHERE (k.alle_freq < .{2:s} OR k.alle_freq IS NULL) LIMIT 5
#             """.format(args.db, gene_statement, args.kav)
#     query
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
#             (SELECT DISTINCT p.sample_id, p.allele1varquality, p.totalreadcount, ens.gene_name, p.chr, p.start, p.stop,
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
#              WHERE (k.alle_freq < .{2:s} OR k.alle_freq IS NULL) LIMIT 5
#             """.format(args.db, gene_statement, args.kav)
#     query
# else:
#     print "Did you select illumina or cgi as your platform? Please check and try again."
#     sys.exit()
#
# #######################
# # connect to database #
# #######################
# print "Connecting to impala..." + "\n"
# from impala.dbapi import connect
#
# conn = connect(host='glados19', port=21050)
# cur = conn.cursor()
#
# #########################
# # Run Query on Impala ##
# #########################
# print "Running query on impala..." + "\n"
# from impala.util import as_pandas
# cur.execute(query)
# query_results = as_pandas(cur)
# if len(query_results) > 0:
#     print "Here's a preview of results found: \n" + (query_results)
# else:
#     print "No results found"

# TESTING #
import pandas as pd
query_df = pd.read_csv('/Users/selasady/impala_scripts/queries/brady_requests/test_results.csv', dtype=str)


############################
# find candidate variants ##
############################
# use variant id to match by gene, chr and position and check inheritance

hom_alt = []
hom_ref = []
comp_het = []

by_variant = query_df.groupby('variant_id')

test = query_df[query_df['gene_name'] == "HSD17B13"]

#print len(test[(test.member == "NB")]) > 0

#     if nb is het
#         and more than one variant per gene
#         and variant from mom and dad are different
#         then mark as comp het


#check for hom_alt and hom_ref rare variants
if (len(test[(test.member == "NB")]) > 0) and (len(test[(test.member == "M")]) > 0) and (len(test[(test.member == "F")]) > 0):
    if (len(test[(test['member']=='NB') & (test['gt']=='0/0')]) > 0 ) \
            and (len(test[(test['member']=='M') & (test['gt']=='0/1')]) > 0 ) \
            and (len(test[(test['member']=='F') & (test['gt']=='0/1')]) > 0 ):
        hom_ref.append(test)
    if (len(test[(test['member']=='NB') & (test['gt']=='1/1')]) > 0 ) \
            and (len(test[(test['member']=='M') & (test['gt']=='0/1')]) > 0 ) \
            and (len(test[(test['member']=='F') & (test['gt']=='0/1')]) > 0 ):
        hom_alt.append(test)

#check for comp hets per gene
if test.variant_id.nunique() > 1:
    print "sweet"



#print test[(test['member']=='NB') & (test['gt']=='0/0')]


# for name, group in by_variant:
#     print name
#     print group





#d.groupby('journey')['mode'].apply(lambda g: 'BUS' in g.values and 'RTS' in g.values)





#
#






#find mie

#
# ########################
# ## Find MIE in Het NB ##
# #########################
# nb_het = nb[which(nb$gt == "0/1"),]
#
# #find matching variants from each parent
# mom_vars = mom[grep(paste(nb_het$variant_id, collapse="|"), mom$variant_id),] #10 vars
# dad_vars = dad[grep(paste(nb_het$variant_id, collapse="|"), dad$variant_id),] #9 vars
#
# #merge together to find intersection of gene:chr:pos
# mie_het_nb_mom = merge(nb_het, mom_vars, by = "variant_id")
# mie_hets_nb_dad = merge(nb_het, dad_vars, by = "variant_id")
# mie_hets = merge(mie_het_nb_mom, mie_hets_nb_dad, by= "variant_id")
#
# #get rid of unnecessary columns
# mie_hets = mie_hets[,c(1:13,25:26,38:39)]
# colnames(mie_hets) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_freqPct",
#                        "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt", "nb_gt",
#                        "m_alt", "m_gt", "f_alt","f_gt")
#
# #variants not in mie_hets and not marked as comp_het are missing info from one parent
# missing_info = nb_het[grep(paste(mie_hets$variant_id, collapse="|"), nb_het$variant_id, invert=TRUE),]
# missing_info = missing_info[grep(paste(comp_het.df$variant_id, collapse="|"), missing_info$variant_id, invert=TRUE),]
#
# #find variants that are in congruence with MI laws
# no_mie = rbind(
#   #if the mother is 0/0 and father is 0/1, then the nb alt must match the father's alt or be null
#   mie_hets[(mie_hets$m_gt == "0/0" & mie_hets$f_gt == "0/1" & (mie_hets$nb_alt == mie_hets$f_alt | mie_hets$nb_alt == "NULL")),],
#   #if the father is 0/0 and the mother is 0/1, then the nb alt must match the mother's alt or be null
#   mie_hets[(mie_hets$f_gt == "0/0" & mie_hets$m_gt == "0/1" & (mie_hets$nb_alt == mie_hets$m_alt| mie_hets$nb_alt == "NULL")),],
#   #if both parents are het, nb alt and parent alts must be the same, or nb alt is null
#   mie_hets[((mie_hets$m_gt == "0/1" & mie_hets$f_gt== "0/1") & (mie_hets$nb_alt == mie_hets$f_alt & mie_hets$nb_alt == mie_hets$m_alt)| (mie_hets$m_alt == mie_hets$f_alt & mie_hets$nb_alt == "NULL")),],
#   #if the mother is 1/1, then the father must be 0/1 or 0/0 and the nb alt must match the mother's alt
#   #or be null
#   mie_hets[(mie_hets$m_gt == "1/1" & (mie_hets$f_gt == "0/1" | mie_hets$f_gt == "0/0") & (mie_hets$nb_alt == mie_hets$m_alt & mie_hets$nb_alt == mie_hets$f_alt | mie_hets$m_alt == mie_hets$f_alt & mie_hets$nb_alt == "NULL")),],
#   #if the father is 1/1, then the mother must be 0/1 or 0/0 and the nb alt must match the fathers's alt
#   mie_hets[(mie_hets$f_gt == "1/1" & (mie_hets$m_gt == "0/1" | mie_hets$m_gt == "0/0") & (mie_hets$nb_alt == mie_hets$f_alt & mie_hets$nb_alt == mie_hets$m_alt | mie_hets$f_alt == mie_hets$m_alt & mie_hets$nb_alt == "NULL")),])
#
# #if mie_het variants are not in no_mie set, then they are potentially MIE
# het_mie = mie_hets[grep(paste(no_mie$variant_id, collapse="|"), mie_hets$variant_id, invert=TRUE),]
#
# ###################################
# ## Find MIE in homozygous ref nb ##
# ###################################
# nb_hom_ref = nb[which(nb$gt == "0/0"),]
#
# #find equivalent variants in parents
# mom_hom_vars = na.omit(mom[match(nb_hom_ref$variant_id, mom$variant_id),])
# dad_hom_vars = na.omit(dad[match(nb_hom_ref$variant_id, dad$variant_id),])
#
# #merge together to find matching parent varients
# mie_hom_nb_mom = merge(nb_hom_ref, mom_hom_vars, by = "variant_id")
# mie_hom_nb_dad = merge(mie_hom_nb_mom, dad_hom_vars, by = "variant_id")
# mie_homs = merge(mie_hom_nb_mom, mie_hom_nb_dad, by= "variant_id")
#
# #clean up results
# mie_homs = mie_homs[,c(1:13,25:26,38:39)]
# colnames(mie_homs) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_freqPct",
#                        "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt", "nb_gt",
#                        "m_alt", "m_gt", "f_alt","f_gt")
#
# #variants that are not in mie_homs and not marked hom ref are missing parent info
# missing_info_miehoms = nb_hom_ref[grep(paste(mie_homs$variant_id, collapse="|"), nb_hom_ref$variant_id, invert=TRUE),]
# missing_info_miehoms = missing_info[grep(paste(hom_ref.df$variant_id, collapse="|"), missing_info$variant_id, invert=TRUE),]
# missing_info_miehoms
#
# #find variants that are in congruence with MI laws
# no_mie_homs = rbind(
#   #if the mother is 0/1, then the dad must be 0/0
#   mie_homs[(mie_homs$m_gt == "0/1" & mie_homs$f_gt == "0/0"),],
#   #if the father is 0/1, then the mother must be 0/0
#   mie_homs[(mie_homs$d_gt == "0/1" & mie_homs$m_gt == "0/0"),],
#   #both parents are 0/0
#   mie_homs[(mie_homs$m_gt == "0/0" & mie_homs$f_gt== "0/0"),]
# )
#
# #if nb_hom variants are not in no_mie_homs set, then they are MIE
# hom_ref_mie = na.omit(nb_hom_ref[!match(nb_hom_ref$variant_id, no_mie_homs$variant_id),])
#
# ###################################
# ## Find MIE in homozygous alt nb ##
# ###################################
# nb_hom_alt = nb[which(nb$gt == "1/1"),]
#
# #find equivalent variants in parents
# mom_hom_alt = na.omit(mom[match(nb_hom_alt$variant_id, mom$variant_id),])
# dad_hom_alt = na.omit(dad[match(nb_hom_alt$variant_id, dad$variant_id),])
#
# #merge together to find intersection of gene:chr:pos
# mie_alt = merge(nb_hom_alt, mom_hom_alt, by = "variant_id")
# mie_alt = merge(mie_alt, dad_hom_alt, by = "variant_id")
#
# #clean up resuts
# mie_alts = mie_alt[,c(1:13,25:26,38:39)]
# colnames(mie_alts) = c("variant_id", "qual", "filter", "rsID", "kav_pct",
#                        "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt",
#                        "nb_gt", "m_alt", "m_gt", "f_alt", "f_gt")
#
# #find variants that are in congruence with MI laws
# mie_hom_alt = rbind(
#   #if the mother is 0/1, then the father must be 1/1
#   mie_alts[(mie_alts$m_gt == "0/1" & mie_alts$f_gt == "1/1"),],
#   #if the father is 0/1, then the mother must be 1/1
#   mie_alts[(mie_alts$d_gt == "0/1" & mie_alts$m_gt == "1/1"),],
#   #both parents are 1/1
#   mie_alts[(mie_alts$m_gt == "1/1" & mie_alts$f_gt == "1/1"),]
# )
#
# #if nb_hom variants are not in no_mie_homs set, then they are MIE
# hom_alt_mie = na.omit(nb_hom_alt[!match(nb_hom_alt$variant_id, mie_hom_alt$variant_id),])
#
# ####################################
# ## Merge Results and save to file ##
# ####################################
# #remove extraneous columns from comp hets
# comp_het_nb.df = comp_het.df[,c(1:13)]
#
# #add parent information
# comp_het_nb.df$m_alt = mom[match(comp_het_nb.df$variant_id, mom$variant_id),]$alt
# comp_het_nb.df$m_gt = mom[match(comp_het_nb.df$variant_id, mom$variant_id),]$gt
# comp_het_nb.df$f_alt = dad[match(comp_het_nb.df$variant_id, dad$variant_id),]$alt
# comp_het_nb.df$f_gt = dad[match(comp_het_nb.df$variant_id, dad$variant_id),]$gt
#
# #make column names the same
# colnames(comp_het_nb.df) = colnames(hom_alt.df)
#
# #add labels to each
# if (dim(hom_alt.df)[1] > 0){
#   hom_alt.df$vartype= "hom_alt"
#   df_list = list("hom_alt.df")
# }
# if (dim(hom_ref.df)[1] > 0){
#   hom_ref.df$vartype = "hom_ref"
#   df_list = c(df_list, "hom_ref.df")
# }
# if (dim(comp_het.df)[1] > 0){
#   comp_het_nb.df$vartype = "comp_het"
#   df_list = c(df_list, "comp_het_nb.df")
# }
# if (dim(het_mie)[1] > 0){
#   het_mie$vartype = "het_mie"
#   df_list = c(df_list, "het_mie")
# }
# if (dim(hom_ref_mie)[1] > 0){
#   hom_ref_mie$vartype = "hom_ref_mie"
#   df_list = c(df_list, "hom_ref_mie")
# }
# if (dim(hom_alt_mie)[1] > 0){
#   hom_alt_mie$vartype = "hom_alt_mie"
#   df_list = c(df_list, "hom_alt_mie")
# }
#
# results.df = unique(do.call("rbind", lapply(df_list, get)))

