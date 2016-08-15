#!/usr/bin/env pyspark


import illumina_snpeff as snp
import sys
import subprocess as sp
from pyspark import SparkContext, SparkConf, SQLContext
from impala.dbapi import connect

###################################
## setup variables and functions ##
###################################

# ITMI options
# ilmn_spark_prefix = 'hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/wgs_ilmn.db/'
ilmn_spark_prefix = 'hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/wgs_ilmn.db/'
ilmn_impala_prefix = 'wgs_ilmn.'
anno_spark_prefix = 'hdfs://ip-10-0-0-118.ec2.internal:8020/itmi/anno_grch37.db/'
anno_impala_prefix = 'anno_grch37.'
appname= 'gv_illumina'
chroms = sorted(map(str, range(1,23) + ["M", "X", "Y"]))
var_blocks = range(0,251)

# configure connection to impala
impala_conn = connect(host="localhost", port=21050, timeout=10000, user='ec2-user')
cur = impala_conn.cursor()

# configure connection to spark
conf = (SparkConf().setAppName(appname))
sc = SparkContext(conf=conf)
sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")
sqlContext.sql("SET spark.sql.parquet.cacheMetadata=true")

# read in table, convert to parquet, store as temp
def impala_query(input_query):
    print("Running query: {}").format(input_query)
    try:
        cur.execute(input_query)
    except Exception as e:
        print (e)

def run_query(input_query):
    sqlContext.sql(input_query)

def check_tables(table1, table2):
    print ("Checking that rows were preserved between {} and {}".format(table1, table2))
    in_table1 = sqlContext.parquetFile("{prefix}/{table}".format(prefix=ilmn_spark_prefix, table=table1))
    in_table1.registerTempTable('t1_df')
    in_table2 = sqlContext.parquetFile("{prefix}/{table}".format(prefix=ilmn_spark_prefix, table=table2))
    in_table2.registerTempTable('t2_df')
    count1 = sqlContext.sql("SELECT COUNT(1) FROM t1_df").collect()
    count2 = sqlContext.sql("SELECT COUNT(1) FROM t2_df").collect()
    if count1 <= count2:
        impala_query("drop table {}{}".format(ilmn_impala_prefix, table1))
    else:
        sys.exit("{} has less rows than {}.".format(table2, table1))

def shut_down():
    sqlContext.clearCache()
    sc.stop()
    cur.close()

########################################################
# create tables needed to store data along the way   ###
########################################################

# create table to store all distinct illumina + ref variants
create_ilmn_vars = "create table {prefix}ilmn_vars \
( \
  var_id string, \
  pos int, \
  ref string, \
  allele string) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix= ilmn_impala_prefix)

create_vars_dbsnp = "  \
create table {prefix}vars_dbsnp \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_vars_kaviar = " \
create table {prefix}vars_kaviar \
( var_id STRING, \
 pos      INT, \
 ref STRING, \
 allele STRING, \
 rs_id STRING, \
 dbsnp_buildid STRING, \
 kav_freq FLOAT, \
 kav_source STRING) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_vars_clinvar = "  \
create table {prefix}vars_clinvar \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string, \
 kav_freq float, \
 kav_source string, \
 clin_sig string, \
 clin_dbn string \
) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_vars_hgmd = "  \
create table {prefix}vars_hgmd \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string, \
 kav_freq float, \
 kav_source string, \
 clin_sig string, \
 clin_dbn string, \
 hgmd_id string, \
 hgmd_varclass string \
) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_vars_dbnsfp = "  \
create table {prefix}vars_dbnsfp \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string, \
 kav_freq float, \
 kav_source string, \
 clin_sig string, \
 clin_dbn string, \
 hgmd_id string, \
 hgmd_varclass string, \
 cadd_raw float, \
 dann_score float, \
 interpro_domain string \
) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_vars_ensembl = " \
create table {prefix}vars_ensembl \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string, \
 kav_freq float, \
 kav_source string, \
 clin_sig string, \
 clin_dbn string, \
 hgmd_id string, \
 hgmd_varclass string, \
 cadd_raw float, \
 dann_score float, \
 interpro_domain string, \
 strand string, \
 gene_name string, \
 gene_id string, \
 transcript_name string, \
 transcript_id string, \
 exon_name string, \
 exon_number int \
) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_vars_coding = " \
create table {prefix}vars_coding \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string, \
 kav_freq float, \
 kav_source string, \
 clin_sig string, \
 clin_dbn string, \
 hgmd_id string, \
 hgmd_varclass string, \
 cadd_raw float, \
 dann_score float, \
 interpro_domain string, \
 strand string, \
 gene_name string, \
 gene_id string, \
 transcript_name string, \
 transcript_id string, \
 exon_name string, \
 exon_number int, \
 effect string, \
 impact string, \
 feature string, \
 feature_id string, \
 biotype string, \
 rank int, \
 hgvs_c string, \
 hgvs_p string \
) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

create_global_vars = " \
create table {prefix}global_vars \
( \
 var_id string, \
 pos int, \
 ref string, \
 allele string, \
 rs_id string, \
 dbsnp_buildid string, \
 kav_freq float, \
 kav_source string, \
 clin_sig string, \
 clin_dbn string, \
 hgmd_id string, \
 hgmd_varclass string, \
 cadd_raw float, \
 dann_score float, \
 interpro_domain string, \
 strand string, \
 gene_name string, \
 gene_id string, \
 transcript_name string, \
 transcript_id string, \
 exon_name string, \
 exon_number int, \
 effect string, \
 impact string, \
 feature string, \
 feature_id string, \
 biotype string, \
 rank int, \
 hgvs_c string, \
 hgvs_p string, \
 var_type string, \
 ppc_rating string \
) \
partitioned by (chrom string, blk_pos int) \
STORED AS PARQUET;".format(prefix=ilmn_impala_prefix)

# create list of tables to create
create_tables_list = [create_ilmn_vars, create_vars_dbsnp, create_vars_kaviar, create_vars_clinvar,
               create_vars_hgmd, create_vars_dbnsfp, create_vars_ensembl, create_vars_coding,
               create_global_vars]

# create each table in the list
# for query in create_tables_list:
#     impala_query(query)

############################################################
# create table of all variants found by joining distinct ###
# illumina variants with distinct reference variants     ###
############################################################

# insert variants into {prefix}ilmn_vars by chromosome and block_pos
for chrom in chroms:
    for pos in var_blocks:
        print ("Inserting chrom {} blk_pos {} into ilmn_vars").format(chrom, pos)
        insert_ilmn_vars = '''
        insert into {ilmn_db}ilmn_vars partition (chrom, blk_pos)
        SELECT var_id, pos, ref, allele, chrom, blk_pos FROM {ilmn_db}vcf_distinct WHERE chrom = '{chrom}' AND blk_pos = {pos}
        UNION
        SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM {anno_db}dbnsfp_distinct_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
        UNION
        SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM {anno_db}kaviar_distinct_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
        UNION
        SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM {anno_db}clinvar_distinct_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
        UNION
        SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM {anno_db}dbsnp_test WHERE chrom = '{chrom}' AND blk_pos = {pos}
        UNION
        SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM {anno_db}hgmd_test WHERE chrom = '{chrom}' AND blk_pos = {pos};
        '''.format(chrom=chrom, pos=pos, ilmn_db=ilmn_spark_prefix, anno_db=anno_spark_prefix)
        # run_query(insert_ilmn_vars)

ilmn_vars_stats = 'compute stats {ilmn_db}ilmn_vars;'.format(ilmn_db=ilmn_impala_prefix)
# impala_query(ilmn_vars_stats)

# ###################
# ### run snpeff  ###
# ###################
# run_snp_cmd = "nohup python illumina_snpeff.py"
# sp.Popen(run_snp_cmd,shell=True)

#################################
### ADD VARIANT  ANNOTATIONS  ###
#################################
# add rsID from dbSNP
for chrom in chroms:
    for pos in var_blocks:
        add_dbsnp = '''
            insert into {ilmn_db}vars_dbsnp partition (chrom, blk_pos)
            with vars as (
            SELECT v.var_id, v.pos, v.ref, v.allele, v.chrom, v.blk_pos
            FROM {ilmn_db}ilmn_vars v
            WHERE v.chrom = '{chrom}'
            AND v.blk_pos = {pos}
            ),
            dbsnp as (
            SELECT d.rs_id, d.dbsnpbuildid as dbsnp_buildid, d.var_id
            from {anno_db}dbsnp d
            )
            SELECT vars.var_id, vars.pos, vars.ref, vars.allele, dbsnp.rs_id, dbsnp.dbsnp_buildid,
                vars.chrom, vars.blk_pos
            from vars
            LEFT JOIN dbsnp
             ON vars.var_id = dbsnp.var_id;
            '''.format(chrom=chrom, pos=pos, ilmn_db=ilmn_spark_prefix, anno_db=anno_spark_prefix)
        # run_query(add_dbsnp)

ilmn_dbsnp_stats = 'compute stats {ilmn_db}vars_dbsnp;'.format(ilmn_db=ilmn_impala_prefix)
# impala_query(ilmn_dbsnp_stats)
# check_tables('ilmn_vars', 'vars_dbsnp')

# add kaviar frequency and source from Kaviar
for chrom in chroms:
    for pos in var_blocks:
        add_kaviar = '''
        insert into {ilmn_db}vars_kaviar partition (chrom, blk_pos)
    WITH vars AS
      (SELECT v.var_id,
             v.pos,
             v.ref,
             v.allele,
             v.rs_id,
             v.dbsnp_buildid,
             v.chrom,
             v.blk_pos
      FROM {ilmn_db}vars_dbsnp v
      WHERE v.chrom = '{chrom}'
      AND v.blk_pos = {pos} ),
      kav as (
        SELECT k.var_id, k.kav_freq, k.kav_source
        FROM {anno_db}kaviar_distinct k
        )
        SELECT vars.var_id, vars.pos, vars.ref, vars.allele, vars.rs_id,
            vars.dbsnp_buildid, kav.kav_freq, kav.kav_source, vars.chrom,
            vars.blk_pos
        FROM vars
        LEFT JOIN kav
         ON vars.var_id = kav.var_id;
        '''.format(chrom=chrom, pos=pos, ilmn_db=ilmn_spark_prefix, anno_db=anno_spark_prefix)
        # run_query(add_kaviar)

# compute stats and check that rows were preserved
ilmn_kaviar_stats = 'compute stats {ilmn_db}vars_kaviar;'.format(ilmn_db=ilmn_impala_prefix)
# impala_query(ilmn_kaviar_stats)
# check_tables('vars_dbsnp', 'vars_kaviar')

# add clinvar significance and disease identification from clinVar
for chrom in chroms:
    for pos in var_blocks:
        add_clinvar = '''
        # insert into {prefix}vars_clinvar partition (chrom, blk_pos)
        with vars as (
                SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
                    v.dbsnp_buildid, v.kav_freq, v.kav_source,
                    v.chrom, v.blk_pos
                FROM {prefix}vars_kaviar v
                WHERE v.chrom = '{chrom}'
                AND v.blk_pos = {pos}
          ),
          clin as (
            SELECT c.var_id, c.clin_sig, c.clin_dbn
            FROM anno_grch37.clinvar_distinct c
            )
            SELECT vars.var_id, vars.pos, vars.ref, vars.allele, vars.rs_id,
                   vars.dbsnp_buildid, vars.kav_freq, vars.kav_source,
                   clin.clin_sig, clin.clin_dbn, vars.chrom, vars.blk_pos
                   FROM vars
        LEFT JOIN clin
         ON vars.var_id = clin.var_id;
        '''.format(chrom=chrom, pos=pos)
        run_query(add_clinvar)

# # compute stats and check that rows were preserved
# run_query("compute stats {prefix}vars_clinvar;")
# check_tables('{prefix}vars_kaviar', '{prefix}vars_clinvar')

# # add hgmd ratings
# for chrom in chroms:
#     for pos in var_blocks:
#         add_hgmd = '''
#         insert into {prefix}vars_hgmd partition (chrom, blk_pos)
#         WITH vars as (
#             SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
#                 v.dbsnp_buildid, v.kav_freq, v.kav_source,
#                 v.clin_sig, v.clin_dbn, v.chrom, v.blk_pos
#             FROM {prefix}vars_clinvar v
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
#         run_query(add_hgmd)
#
# # compute stats
# run_query("compute stats {prefix}vars_hgmd;")
#
# check_tables('{prefix}vars_clinvar', '{prefix}vars_hgmd')

# # add cadd, dann and interpro domain from dbnsfp
# for chrom in chroms:
#     for pos in var_blocks:
#         add_dbnsfp = '''
#         insert into {prefix}vars_dbnsfp partition (chrom, blk_pos)
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
#           FROM {prefix}vars_hgmd v
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
#         run_query(add_dbnsfp)
#
# # compute stats
# run_query("compute stats {prefix}vars_dbnsfp;")
#
# check_tables('{prefix}vars_hgmd', '{prefix}vars_dbnsfp')
#
#
# # add gene, transcript and exon id and names from ensembl
# for chrom in chroms:
#     for pos in var_blocks:
#         add_ensembl = '''
#         INSERT INTO TABLE {prefix}vars_ensembl partition(chrom, blk_pos)
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
#           FROM {prefix}vars_dbnsfp v
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
#         run_query(add_ensembl)
#
# check_tables('{prefix}vars_dbnsfp', '{prefix}vars_ensembl')
# run_query("compute stats {prefix}vars_ensembl")
#
# #######################################
# ### ADD snpEff Coding Consequences  ###
# #######################################
# # create table to store snpeff results
# create_snpeff_table = '''
# create external table {prefix}snpeff_results
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
# '''.format(hdfs_out)
#
# # make sure that the snpeff pipeline has finished running
# snpeff_thread.join()
#
# # add snpeff annotations
# run_query(create_snpeff_table)
# run_query("compute stats {prefix}snpeff_results")
#
# # create partitioned table
# create_snpeff_partitioned = '''
# create table {prefix}snpeff_partitioned
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
# run_query(create_snpeff_partitioned)
#
# # insert results into partitioned table
# for chrom in chroms:
#     for pos in var_blocks:
#         insert_snpeff_partitioned = '''
#         insert into table {prefix}snpeff_partitioned  partition (chrom, blk_pos)
#         select var_id, pos, ref, alt,gene,gene_id,affect,
#         impact, feature, feature_id, biotype, rank,
#         distance, hgvs_c, hgvs_p, chrom, blk_pos
#         from {prefix}snpeff_results
#         where chrom = '{chrom}'
#         and blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         run_query(insert_snpeff_partitioned)
#
# run_query("compute stats {prefix}snpeff_partitioned")
#
# # join snpeff results with annotated variants
# for chrom in chroms:
#     for pos in var_blocks:
#         add_coding = '''
#         insert into table {prefix}vars_coding partition(chrom, blk_pos)
#         select v.var_id,v.pos,v.ref,v.allele,v.rs_id,
#              v.dbsnp_buildid,v.kav_freq,v.kav_source,
#              v.clin_sig,v.clin_dbn,v.hgmd_id, v.hgmd_varclass,
#              v.cadd_raw, v.dann_score,v.interpro_domain,
#              v.strand,v.gene_name,v.gene_id,v.transcript_name,
#              v.transcript_id,v.exon_name,v.exon_number,
#              s.affect,s.impact,s.feature,s.feature_id,s.biotype,
#              s.rank,s.hgvs_c,s.hgvs_p, v.chrom, v.blk_pos
#          FROM {prefix}vars_ensembl v
#          LEFT JOIN {prefix}snpeff_partitioned s
#                        ON v.var_id = s.var_id
#          WHERE v.chrom = '{chrom}' and s.chrom = '{chrom}'
#          AND v.blk_pos = {pos} and s.blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         run_query(add_coding)
#
# run_query("compute stats {prefix}vars_coding")
#
# ####################################
# #### ANNOTATE WITH PPC NOTATION  ###
# ####################################
#
# # add ppc notation to create final global var table
# create_gv = '''
# create table {prefix}global_vars
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
# run_query(create_gv)
#
# for chrom in chroms:
#     for pos in var_blocks:
#         add_ppc = '''
#         insert into {prefix}global_vars partition (chrom, blk_pos)
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
#         FROM {prefix}vars_coding
#         WHERE chrom = '{chrom}'
#         AND blk_pos = {pos};
#         '''.format(chrom=chrom, pos=pos)
#         run_query(add_ppc)
#
# run_query("compute stats {prefix}global_vars")
#
# print("Global variants pipeline complete.")

# clean up temp and close connection
shut_down()