#!/usr/bin/env python
import illumina_snpeff as snp
import sys
import subprocess as sp

# ITMI options
impala_host = 'localhost'
impala_port = 21050
impala_user_name = 'ec2-user'
hdfs_path = 'elasasu/'
vcf_dir = '/home/ec2-user/elasasu/impala_scripts/global_vars/illumina_gv'

# chrom_list = map(str, range(1, 23)) + ['X', 'Y', 'M']
chrom_list = ['M']

# instantiate snpeff script and variables
snpeff = snp.snpeff(vcf_dir, impala_host, impala_port, impala_user_name, hdfs_path)

########################################################
# create tables needed to store data along the way   ###
########################################################

# create table to store all distinct illumina + ref variants
create_ilmn_vars = '''
create table wgs_ilmn.ilmn_vars
(
  var_id string,
  pos int,
  ref string,
  allele string)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

create_vars_dbsnp = '''
create table wgs_ilmn.vars_dbsnp
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

create_vars_kaviar = '''
create table wgs_ilmn.vars_kaviar
(
 var_id STRING,
 pos      INT,
 ref STRING,
 allele STRING,
 rs_id STRING,
 dbsnp_buildid STRING,
 kav_freq FLOAT,
 kav_source STRING
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

create_vars_clinvar = '''
create table wgs_ilmn.vars_clinvar
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string,
 kav_freq float,
 kav_source string,
 clin_sig string,
 clin_dbn string
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''
create_vars_hgmd = '''
create table wgs_ilmn.vars_hgmd
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string,
 kav_freq float,
 kav_source string,
 clin_sig string,
 clin_dbn string,
 hgmd_id string,
 hgmd_varclass string
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

create_vars_dbnsfp = '''
create table wgs_ilmn.vars_dbnsfp
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string,
 kav_freq float,
 kav_source string,
 clin_sig string,
 clin_dbn string,
 hgmd_id string,
 hgmd_varclass string,
 cadd_raw float,
 dann_score float,
 interpro_domain string
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''
create_vars_ensembl = '''
create table wgs_ilmn.vars_ensembl
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string,
 kav_freq float,
 kav_source string,
 clin_sig string,
 clin_dbn string,
 hgmd_id string,
 hgmd_varclass string,
 cadd_raw float,
 dann_score float,
 interpro_domain string,
 strand string,
 gene_name string,
 gene_id string,
 transcript_name string,
 transcript_id string,
 exon_name string,
 exon_number int
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

create_vars_coding = '''
create table wgs_ilmn.vars_coding
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string,
 kav_freq float,
 kav_source string,
 clin_sig string,
 clin_dbn string,
 hgmd_id string,
 hgmd_varclass string,
 cadd_raw float,
 dann_score float,
 interpro_domain string,
 strand string,
 gene_name string,
 gene_id string,
 transcript_name string,
 transcript_id string,
 exon_name string,
 exon_number int,
 effect string,
 impact string,
 feature string,
 feature_id string,
 biotype string,
 rank int,
 hgvs_c string,
 hgvs_p string
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

create_global_vars = '''
create table wgs_ilmn.global_vars
(
 var_id string,
 pos int,
 ref string,
 allele string,
 rs_id string,
 dbsnp_buildid string,
 kav_freq float,
 kav_source string,
 clin_sig string,
 clin_dbn string,
 hgmd_id string,
 hgmd_varclass string,
 cadd_raw float,
 dann_score float,
 interpro_domain string,
 strand string,
 gene_name string,
 gene_id string,
 transcript_name string,
 transcript_id string,
 exon_name string,
 exon_number int,
 effect string,
 impact string,
 feature string,
 feature_id string,
 biotype string,
 rank int,
 hgvs_c string,
 hgvs_p string,
 var_type string,
 ppc_rating string
)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;
'''

# create list of tables to create
create_tables_list = [create_vars_dbsnp, create_vars_kaviar, create_vars_clinvar,
               create_vars_hgmd, create_vars_dbnsfp, create_vars_ensembl, create_vars_coding,
               create_global_vars]

# create each table in the list
# for query in create_tables_list:
#     snpeff.run_query(query)

############################################################
# create table of all variants found by joining distinct ###
# illumina variants with distinct reference variants     ###
############################################################

# insert variants into wgs_ilmn.ilmn_vars by chromosome and block_pos
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         print ("Running query for chrom {} blk_pos {}").format(chrom, pos)
#         insert_ilmn_vars = '''
#         insert into wgs_ilmn.ilmn_vars partition (chrom, blk_pos)
#         SELECT var_id, pos, ref, allele, chrom, blk_pos FROM wgs_ilmn.test_vars WHERE chrom = '{chrom}' AND blk_pos = {pos}
#         UNION
#         SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.dbnsfp_distinct_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
#         UNION
#         SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.kaviar_distinct_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
#         UNION
#         SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.clinvar_distinct_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
#         UNION
#         SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.dbsnp_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
#         UNION
#         SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.hgmd_test WHERE chrom = '{chrom}' AND blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         # snpeff.run_query(insert_ilmn_vars)
#
# # snpeff.run_query("compute stats wgs_ilmn.ilmn_vars;")
#
# ###################
# ### run snpeff  ###
# ###################
# run_snp_cmd = "nohup python illumina_snpeff.py"
# sp.Popen(run_snp_cmd,shell=True)

#################################
### ADD VARIANT  ANNOTATIONS  ###
#################################

# check that each table has at least as many rows as the previous table
# if new table passes check, drop previous table to save disk space
def check_tables(table1, table2):
    count1 = snpeff.pandas_query("SELECT COUNT(1) FROM {}".format(table1))
    count2 = snpeff.pandas_query("SELECT COUNT(1) FROM {}".format(table2))
    if int(count2.ix[0]) >= int(count1.ix[0]) >= 1000:
        snpeff.run_query("drop table {}".format(table1))
    else:
        sys.exit("{} has less rows than {}.".format(table2, table1))

# add rsID from dbSNP
for chrom in chrom_list:
    for pos in snpeff.var_blocks:
        add_dbsnp = '''
            insert into wgs_ilmn.vars_dbsnp partition (chrom, blk_pos)
            with vars as (
            SELECT v.var_id, v.pos, v.ref, v.allele, v.chrom, v.blk_pos
            FROM wgs_ilmn.ilmn_vars v
            WHERE v.chrom = '{chrom}'
            AND v.blk_pos = {pos}
            ),
            dbsnp as (
            SELECT d.rs_id, d.dbsnpbuildid as dbsnp_buildid, d.var_id
            from anno_grch37.dbsnp d
            )
            SELECT vars.var_id, vars.pos, vars.ref, vars.allele, dbsnp.rs_id, dbsnp.dbsnp_buildid,
                vars.chrom, vars.blk_pos
            from vars
            LEFT JOIN dbsnp
             ON vars.var_id = dbsnp.var_id;
            '''.format(chrom=chrom, pos=pos)
        # snpeff.run_query(add_dbsnp)

# compute stats
# snpeff.run_query("compute stats wgs_ilmn.vars_dbsnp;")

# add kaviar frequency and source from Kaviar
for chrom in chrom_list:
    for pos in snpeff.var_blocks:
        add_kaviar = '''
        insert into wgs_ilmn.vars_kaviar partition (chrom, blk_pos)
    WITH vars AS
      (SELECT v.var_id,
             v.pos,
             v.ref,
             v.allele,
             v.rs_id,
             v.dbsnp_buildid,
             v.chrom,
             v.blk_pos
      FROM wgs_ilmn.vars_dbsnp v
      WHERE v.chrom = '{chrom}'
      AND v.blk_pos = {pos} ),
      kav as (
        SELECT k.var_id, k.kav_freq, k.kav_source
        FROM anno_grch37.kaviar_distinct k
        )
        SELECT vars.var_id, vars.pos, vars.ref, vars.allele, vars.rs_id,
            vars.dbsnp_buildid, kav.kav_freq, kav.kav_source, vars.chrom,
            vars.blk_pos
        FROM vars
        LEFT JOIN kav
         ON vars.var_id = kav.var_id;
        '''.format(chrom=chrom, pos=pos)
        snpeff.run_query(add_kaviar)

# compute stats
snpeff.run_query("compute stats wgs_ilmn.vars_kaviar;")

# check_tables('wgs_ilmn.vars_dbsnp', 'wgs_ilmn.vars_kaviar')

# # add clinvar significance and disease identification from clinVar
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         add_clinvar = '''
#         insert into wgs_ilmn.vars_clinvar partition (chrom, blk_pos)
#         with vars as (
#                 SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
#                     v.dbsnp_buildid, v.kav_freq, v.kav_source,
#                     v.chrom, v.blk_pos
#                 FROM wgs_ilmn.vars_kaviar v
#                 WHERE v.chrom = '{chrom}'
#                 AND v.blk_pos = {pos}
#           ),
#           clin as (
#             SELECT c.var_id, c.clin_sig, c.clin_dbn
#             FROM anno_grch37.clinvar_distinct c
#             )
#             SELECT vars.var_id, vars.pos, vars.ref, vars.allele, vars.rs_id,
#                    vars.dbsnp_buildid, vars.kav_freq, vars.kav_source,
#                    clin.clin_sig, clin.clin_dbn, vars.chrom, vars.blk_pos
#                    FROM vars
#         LEFT JOIN clin
#          ON vars.var_id = clin.var_id;
#         '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(add_clinvar)
#
# # compute stats
# snpeff.run_query("compute stats wgs_ilmn.vars_clinvar;")
#
# check_tables('wgs_ilmn.vars_kaviar', 'wgs_ilmn.vars_clinvar')
#
# # add hgmd ratings
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         add_hgmd = '''
#         insert into wgs_ilmn.vars_hgmd partition (chrom, blk_pos)
#         WITH vars as (
#             SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
#                 v.dbsnp_buildid, v.kav_freq, v.kav_source,
#                 v.clin_sig, v.clin_dbn, v.chrom, v.blk_pos
#             FROM wgs_ilmn.vars_clinvar v
#             WHERE v.chrom = '{chrom}'
#             AND v.blk_pos = {pos}
#           ), hgmd as (
#           SELECT h.var_id, h.id as hgmd_id,
#                 h.var_class as hgmd_varclass
#             FROM anno_grch37.hgmd h
#             )
#          SELECT vars.var_id, vars.pos, vars.ref, vars.allele, vars.rs_id,
#                 vars.dbsnp_buildid, vars.kav_freq, vars.kav_source,
#                 vars.clin_sig, vars.clin_dbn, hgmd.hgmd_id, hgmd.hgmd_varclass,
#                 vars.chrom, vars.blk_pos
#             FROM vars
#             LEFT JOIN hgmd
#                ON vars.var_id = hgmd.var_id;
#         '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(add_hgmd)
#
# # compute stats
# snpeff.run_query("compute stats wgs_ilmn.vars_hgmd;")
#
# check_tables('wgs_ilmn.vars_clinvar', 'wgs_ilmn.vars_hgmd')
#
#
# # add cadd, dann and interpro domain from dbnsfp
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         add_dbnsfp = '''
#         insert into wgs_ilmn.vars_dbnsfp partition (chrom, blk_pos)
#         WITH vars AS
#           (SELECT v.var_id,
#                  v.pos,
#                  v.ref,
#                  v.allele,
#                  v.rs_id,
#                  v.dbsnp_buildid,
#                  v.kav_freq,
#                  v.kav_source,
#                  v.clin_sig,
#                  v.clin_dbn,
#                  v.hgmd_id,
#                  v.hgmd_varclass,
#                  v.chrom,
#                  v.blk_pos
#           FROM wgs_ilmn.vars_hgmd v
#           WHERE v.chrom = '{chrom}'
#                   AND v.blk_pos = {pos} ), dbnsfp AS
#           (SELECT d.var_id,
#                  d.cadd_raw,
#                  d.dann_score,
#                  d.interpro_domain
#           FROM anno_grch37.dbnsfp_distinct d )
#         SELECT vars.var_id,
#                  vars.pos,
#                  vars.ref,
#                  vars.allele,
#                  vars.rs_id,
#                  vars.dbsnp_buildid,
#                  vars.kav_freq,
#                  vars.kav_source,
#                  vars.clin_sig,
#                  vars.clin_dbn,
#                  vars.hgmd_id,
#                  vars.hgmd_varclass,
#                  dbnsfp.cadd_raw,
#                  dbnsfp.dann_score,
#                  dbnsfp.interpro_domain,
#                  vars.chrom,
#                  vars.blk_pos
#         FROM vars
#         LEFT JOIN dbnsfp
#             ON vars.var_id = dbnsfp.var_id;
#         '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(add_dbnsfp)
#
# # compute stats
# snpeff.run_query("compute stats wgs_ilmn.vars_dbnsfp;")
#
# check_tables('wgs_ilmn.vars_hgmd', 'wgs_ilmn.vars_dbnsfp')
#
#
# # add gene, transcript and exon id and names from ensembl
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         add_ensembl = '''
#         INSERT INTO TABLE wgs_ilmn.vars_ensembl partition(chrom, blk_pos)
#         WITH vars AS
#           (SELECT v.var_id,
#                  v.pos,
#                  v.ref,
#                  v.allele,
#                  v.rs_id,
#                  v.dbsnp_buildid,
#                  v.kav_freq,
#                  v.kav_source,
#                  v.clin_sig,
#                  v.clin_dbn,
#                  v.hgmd_id,
#                  v.hgmd_varclass,
#                  v.cadd_raw,
#                  v.dann_score,
#                  v.interpro_domain,
#                  v.chrom,
#                  v.blk_pos
#           FROM wgs_ilmn.vars_dbnsfp v
#           WHERE v.chrom = '{chrom}'
#                   AND v.blk_pos = {pos} ), ens AS
#           (SELECT e.strand,
#                  e.gene_name,
#                  e.gene_id,
#                  e.transcript_name,
#                  e.transcript_id,
#                  e.exon_name,
#                  e.exon_number,
#                  e.chrom,
#                  e.pos,
#                  e.stop
#           FROM anno_grch37.ensembl_distinct e )
#         SELECT vars.var_id,
#                  vars.pos,
#                  vars.ref,
#                  vars.allele,
#                  vars.rs_id,
#                  vars.dbsnp_buildid,
#                  vars.kav_freq,
#                  vars.kav_source,
#                  vars.clin_sig,
#                  vars.clin_dbn,
#                  vars.hgmd_id,
#                  vars.hgmd_varclass,
#                  vars.cadd_raw,
#                  vars.dann_score,
#                  vars.interpro_domain,
#                  ens.strand,
#                  ens.gene_name,
#                  ens.gene_id,
#                  ens.transcript_name,
#                  ens.transcript_id,
#                  ens.exon_name,
#                  ens.exon_number,
#                  vars.chrom,
#                  vars.blk_pos
#         FROM vars
#         LEFT JOIN ens
#             ON vars.chrom = ens.chrom
#                 AND (vars.pos
#             BETWEEN ens.pos
#                 AND ens.stop);
#
#
#             '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(add_ensembl)
#
# check_tables('wgs_ilmn.vars_dbnsfp', 'wgs_ilmn.vars_ensembl')
# snpeff.run_query("compute stats wgs_ilmn.vars_ensembl")
#
# #######################################
# ### ADD snpEff Coding Consequences  ###
# #######################################
# # create table to store snpeff results
# create_snpeff_table = '''
# create external table wgs_ilmn.snpeff_results
# (
#   pos int,
#   var_id string,
#   ref string,
#   alt string,
#   gene string,
#   gene_id string,
#   affect string,
#   impact string,
#   feature string,
#   feature_id string,
#   biotype string,
#   rank int,
#   distance string,
#   hgvs_c string,
#   hgvs_p string,
#   chrom string,
#   blk_pos int
#   )
# row format delimited fields terminated by '\t'
# location '{}';
# '''.format(snpeff.hdfs_out)
#
# # make sure that the snpeff pipeline has finished running
# snpeff_thread.join()
#
# # add snpeff annotations
# snpeff.run_query(create_snpeff_table)
# snpeff.run_query("compute stats wgs_ilmn.snpeff_results")
#
# # create partitioned table
# create_snpeff_partitioned = '''
# create table wgs_ilmn.snpeff_partitioned
# (
#   var_id string,
#   pos int,
#   ref string,
#   alt string,
#   gene string,
#   gene_id string,
#   affect string,
#   impact string,
#   feature string,
#   feature_id string,
#   biotype string,
#   rank int,
#   distance string,
#   hgvs_c string,
#   hgvs_p string
#   )
# PARTITIONED BY (chrom string, blk_pos int)
# STORED AS PARQUET;
# '''
#
# snpeff.run_query(create_snpeff_partitioned)
#
# # insert results into partitioned table
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         insert_snpeff_partitioned = '''
#         insert into table wgs_ilmn.snpeff_partitioned  partition (chrom, blk_pos)
#         select var_id, pos, ref, alt,gene,gene_id,affect,
#         impact, feature, feature_id, biotype, rank,
#         distance, hgvs_c, hgvs_p, chrom, blk_pos
#         from wgs_ilmn.snpeff_results
#         where chrom = '{chrom}'
#         and blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(insert_snpeff_partitioned)
#
# snpeff.run_query("compute stats wgs_ilmn.snpeff_partitioned")
#
# # join snpeff results with annotated variants
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         add_coding = '''
#         insert into table wgs_ilmn.vars_coding partition(chrom, blk_pos)
#         select v.var_id,v.pos,v.ref,v.allele,v.rs_id,
#              v.dbsnp_buildid,v.kav_freq,v.kav_source,
#              v.clin_sig,v.clin_dbn,v.hgmd_id, v.hgmd_varclass,
#              v.cadd_raw, v.dann_score,v.interpro_domain,
#              v.strand,v.gene_name,v.gene_id,v.transcript_name,
#              v.transcript_id,v.exon_name,v.exon_number,
#              s.affect,s.impact,s.feature,s.feature_id,s.biotype,
#              s.rank,s.hgvs_c,s.hgvs_p, v.chrom, v.blk_pos
#          FROM wgs_ilmn.vars_ensembl v
#          LEFT JOIN wgs_ilmn.snpeff_partitioned s
#                        ON v.var_id = s.var_id
#          WHERE v.chrom = '{chrom}' and s.chrom = '{chrom}'
#          AND v.blk_pos = {pos} and s.blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(add_coding)
#
# snpeff.run_query("compute stats wgs_ilmn.vars_coding")
#
# ####################################
# #### ANNOTATE WITH PPC NOTATION  ###
# ####################################
#
# # add ppc notation to create final global var table
# create_gv = '''
# create table wgs_ilmn.global_vars
# (
#  var_id string,
#  pos int,
#  ref string,
#  allele string,
#  rs_id string,
#  dbsnp_buildid string,
#  kav_freq float,
#  kav_source string,
#  clin_sig string,
#  clin_dbn string,
#  hgmd_id string,
#  hgmd_varclass string,
#  cadd_raw float,
#  dann_score float,
#  interpro_domain string,
#  strand string,
#  gene_name string,
#  gene_id string,
#  transcript_name string,
#  transcript_id string,
#  exon_name string,
#  exon_number int,
#  effect string,
#  impact string,
#  feature string,
#  feature_id string,
#  biotype string,
#  rank int,
#  hgvs_c string,
#  hgvs_p string,
#  var_type string,
#  ppc_rating string
# )
# partitioned by (chrom string, blk_pos int)
# STORED AS PARQUET;
# '''
#
# snpeff.run_query(create_gv)
#
# for chrom in chrom_list:
#     for pos in snpeff.var_blocks:
#         add_ppc = '''
#         insert into wgs_ilmn.global_vars partition (chrom, blk_pos)
#         SELECT var_id,pos,ref,allele,rs_id,dbsnp_buildid,kav_freq,
#         kav_source, clin_sig, clin_dbn, hgmd_id, hgmd_varclass,
#         cadd_raw, dann_score, interpro_domain, strand,
#         gene_name, gene_id, transcript_name, transcript_id,
#         exon_name, exon_number, effect, impact, feature,
#         feature_id, biotype, rank, hgvs_c, hgvs_p,
#            CASE
#            WHEN (length(ref) > length(allele)) THEN 'deletion'
#            WHEN (length(ref) < length(allele)) THEN 'insertion'
#            WHEN (length(ref) = length(allele)) THEN 'snv'
#            ELSE 'other' END AS 'var_type',
#            CASE
#            WHEN ((length(ref) = length(allele)) and (effect IN ('codings_sequence_variant', 'chromosome',
#            'missense_variant', 'initator_codon_variant', 'stop_retained_variant', 'transcript_variant'))) THEN 'non-synonymous'
#            WHEN ((length(ref) > length(allele)) and (effect = 'frameshift_variant')) THEN 'frameshift-non-synonymous'
#            WHEN ((length(ref) > length(allele)) and (effect = 'disruptive_inframe_insertion')) THEN 'non-synonymous'
#            WHEN ((length(ref) < length(allele)) and (effect = 'frameshift_variant')) THEN 'frameshift-non-synonymous'
#            WHEN ((length(ref) < length(allele)) and (effect = 'disruptive_inframe_deletion')) THEN 'non-synonymous'
#            ELSE 'synonymous' END AS ppc_rating, chrom, blk_pos
#         FROM wgs_ilmn.vars_coding
#         WHERE chrom = '{chrom}'
#         AND blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         snpeff.run_query(add_ppc)
#
# snpeff.run_query("compute stats wgs_ilmn.global_vars")
#
# print("Global variants pipeline complete.")
