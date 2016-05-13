#!/usr/bin/env python
import snpeff as snp

# user args
vcf_dir = '/titan/ITMI1/workspaces/users/selasady/impala_scripts/annotation/snpeff'
impala_host = 'glados14'
impala_port = 21050
impala_user_name = 'selasady'
hdfs_path = '/user/selasady/'

# instantiate snpeff script and variables
snpeff = snp.run_snpeff(vcf_dir, impala_host, impala_port, impala_user_name, hdfs_path)


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
create_tables_list = [create_ilmn_vars, create_vars_dbsnp, create_vars_kaviar, create_vars_clinvar,
               create_vars_hgmd, create_vars_dbnsfp, create_vars_ensembl, create_vars_coding,
               create_global_vars]

# create each table in the list
for query in create_tables_list:
    snpeff.run_query(query)

############################################################
# create table of all variants found by joining distinct ###
# illumina variants with distinct reference variants     ###
############################################################

# insert variants into wgs_ilmn.ilmn_vars by chromosome and block_pos
# for chrom in snpeff.chroms:
#     for pos in snpeff.blk_poss:
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
#         snpeff.run_query(insert_ilmn_vars)

# snpeff.run_query("compute stats wgs_ilmn.ilmn_vars;")

###################
### run snpeff  ###
###################
# snpeff.run_pipeline()

#################################
### ADD VARIANT  ANNOTATIONS  ###
#################################

for chrom in snpeff.chroms:
    for pos in snpeff.blk_pos:
        add_dbsnp = '''
            insert into wgs_ilmn.vars_dbsnp partition (chrom, blk_pos)
            SELECT v.var_id, v.pos, v.ref, v.allele, d.rs_id,
            d.dbsnpbuildid as dbsnp_buildid, v.chrom, v.blk_pos
            FROM wgs_ilmn.ilmn_vars v
            LEFT JOIN anno_grch37.dbsnp d
             ON v.var_id = d.var_id
            WHERE v.chrom = '{chrom}' and d.chrom = '{chrom}'
            AND v.blk_pos = {pos} and d.blk_pos = {pos};
            '''.format(chrom=chrom, pos=pos)

        add_kaviar = '''
        insert into wgs_ilmn.vars_kaviar partition (chrom, blk_pos)
        SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
            v.dbsnp_buildid, k.kav_freq, k.kav_source,
            v.chrom, v.blk_pos
        FROM wgs_ilmn.vars_dbsnp v
        LEFT JOIN anno_grch37.kaviar_distinct k
         ON v.var_id = k.var_id
        WHERE v.chrom = '{chrom}' and k.chrom = '{chrom}'
        AND v.blk_pos = {pos} and k.blk_pos = {pos};
        '''.format(chrom=chrom, pos=pos)

        add_clinvar = '''
        insert into wgs_ilmn.vars_clinvar partition (chrom, blk_pos)
        SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
            v.dbsnp_buildid, v.kav_freq, v.kav_source,
            c.clin_sig, c.clin_dbn, v.chrom, v.blk_pos
        FROM wgs_ilmn.vars_kaviar v
        LEFT JOIN anno_grch37.clinvar_distinct c
         ON v.var_id = c.var_id
        WHERE v.chrom = '{chrom}' and c.chrom = '{chrom}'
        AND v.blk_pos = {pos} and c.blk_pos = {pos};
        '''.format(chrom=chrom, pos=pos)

        add_hgmd = '''
        insert into wgs_ilmn.vars_hgmd partition (chrom, blk_pos)
        SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
            v.dbsnp_buildid, v.kav_freq, v.kav_source,
            v.clin_sig, v.clin_dbn, h.id as hgmd_id,
            h.var_class as hgmd_varclass, v.chrom, v.blk_pos
        FROM wgs_ilmn.vars_clinvar v
        LEFT JOIN anno_grch37.hgmd h
           ON v.var_id = h.var_id
        WHERE v.chrom = '{chrom}' and h.chrom = '{chrom}'
        AND v.blk_pos = {pos} and h.blk_pos = {pos};
        '''.format(chrom=chrom, pos=pos)

        add_dbnsfp = '''
        insert into wgs_ilmn.vars_dbnsfp partition (chrom, blk_pos)
        SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
            v.dbsnp_buildid, v.kav_freq, v.kav_source,
            v.clin_sig, v.clin_dbn, v.hgmd_id,
            v.hgmd_varclass, d.cadd_raw, d.dann_score, d.interpro_domain,
            v.chrom, v.blk_pos
        FROM wgs_ilmn.vars_hgmd v
        LEFT JOIN anno_grch37.dbnsfp_distinct d
        ON v.var_id = d.var_id
        WHERE v.chrom = '{chrom}' and d.chrom = '{chrom}'
        AND v.blk_pos = {pos} and d.blk_pos = {pos};
        '''.format(chrom=chrom, pos=pos)

        add_ensembl = '''
            SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
                v.dbsnp_buildid, v.kav_freq, v.kav_source,
                v.clin_sig, v.clin_dbn, v.hgmd_id,
                v.hgmd_varclass, v.cadd_raw, v.dann_score, v.interpro_domain,
                e.strand, e.gene_name, e.gene_id, e.transcript_name, e.transcript_id,
                e.exon_name, e.exon_number, v.chrom, v.blk_pos
            FROM wgs_ilmn.vars_dbsnfp v
            LEFT JOIN anno_grch37.ensembl_distinct e
               ON v.chrom = e.chrom
               AND (v.pos BETWEEN e.pos and e.stop)
            WHERE v.chrom = '{chrom}' and e.chrom = '{chrom}'
            AND v.blk_pos = {pos} and e.blk_pos = {pos};
            '''.format(chrom=chrom, pos=pos)

        annot_list = [add_dbsnp, add_kaviar, add_clinvar, add_hgmd, add_dbnsfp, add_ensembl]

        for query in annot_list:
            print ("Running query: \n {} \n").format(query)
            snpeff.run_query(query)

# compute stats
snpeff.run_query("compute stats wgs_ilmn.vars_dbnsfp;")

# clean up temp tables
drop_list1 = ['wgs_ilmn.vars_dbsnp', 'wgs_ilmn.vars_kaviar', 'wgs_ilmn.vars_clinvar',
              'wgs_ilmn.vars_hgmd']

for query in drop_list1:
    snpeff.run_query("drop table {}".format(query))

#######################################
### ADD snpEff Coding Consequences  ###
#######################################
# create table to store snpeff results
create_snpeff_table = '''
create external table wgs_ilmn.snpeff_results
(
  chrom string,
  pos int,
  id string,
  ref string,
  alt string,
  gene string,
  gene_id string,
  affect string,
  impact string,
  feature string,
  feature_id string,
  biotype string,
  rank string,
  distance string,
  hgvs_c string,
  hgvs_p string
  )
location '{}/snpeff';
'''.format(hdfs_path)

snpeff.run_query(create_snpeff_table)

# create partitioned table
create_snpeff_partitioned = '''
create table wgs_ilmn.snpeff_partitioned
(
  var_id string,
  pos int,
  id string,
  ref string,
  alt string,
  gene string,
  gene_id string,
  affect string,
  impact string,
  feature string,
  feature_id string,
  biotype string,
  rank string,
  distance string,
  hgvs_c string,
  hgvs_p string
  )
PARTITIONED BY (chrom string, blk_pos int)
STORED AS PARQUET;
'''.format(hdfs_path)

snpeff.run_query(create_snpeff_partitioned)

# insert results into partitioned table
for chrom in snpeff.chroms:
    for pos in snpeff.blk_pos:
        insert_snpeff_partitioned = '''
        insert into table wgs_ilmn.snpeff_partitioned  partition (chrom, blk_pos)
        select concat(chrom, ':', cast(pos as string), ':', ref, alt) as var_id, pos, id, ref, alt,gene,gene_id,affect,
        impact, feature, feature_id, biotype, rank,
        distance, hgvs_c, hgvs_p, chrom, cast(pos/1000000 as int) as blk_pos
        from wgs_ilmn.snpeff_results;
        compute stats wgs_ilmn.snpeff_partitioned
        where chrom = '{chrom}'
        and blk_pos = {pos};
        '''.format(chrom, pos)
        snpeff.run_query(insert_snpeff_partitioned)

snpeff.run_query("compute stats wgs_ilmn.snpeff_partitioned")

# join snpeff results with annotated variants
for chrom in snpeff.chroms:
    for pos in snpeff.blk_pos:
        add_coding = '''
        insert into table wgs_ilmn.vars_coding
        select v.var_id,v.pos,v.ref,v.allele,v.rs_id,
             v.dbsnp_buildid,v.kav_freq,v.kav_source,
             v.clin_sig,v.clin_dbn,v.hgmd_id, v.hgmd_varclass,
             v.cadd_raw, v.dann_score,v.interpro_domain,
             v.strand,v.gene_name,v.gene_id,v.transcript_name,
             v.transcript_id,v.exon_name,v.exon_number,
             s.effect,s.impact,s.feature,s.feature_id,s.biotype,
             s.rank,s.hgvs_c,s.hgvs_p, v.chrom, v.blk_pos
         FROM wgs_ilmn.vars_dbnsfp v
         LEFT JOIN wgs_ilmn.snpeff_partitioned s
                       ON v.var_id = s.var_id
         WHERE v.chrom = '{chrom}' and e.chrom = '{chrom}'
         AND v.blk_pos = {pos} and e.blk_pos = {pos};
         compute stats wgs_ilmn.vars_coding;
        '''.format(chrom, pos)
        snpeff.run_query(add_coding)

####################################
#### ANNOTATE WITH PPC NOTATION  ###
####################################

# add ppc notation to create final global var table
create_gv = '''
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
partitioned by (chrom string, pos_block int)
STORED AS PARQUET;
'''

snpeff.run_query(create_gv)

for chrom in snpeff.chroms:
    for pos in snpeff.blk_pos:
        add_ppc = '''
        insert into wgs_ilmn.global_vars partition (chrom, pos_block)
        SELECT v.*,
           CASE WHEN ((length(ref) == length(allele)) and (effect IN ('codings_sequence_variant', 'chromosome',
           'missense_variant', 'initator_codon_variant', 'stop_retained_variant', 'transcript_variant')) THEN 'non-synonymous',
           WHEN ((length(ref) > length(allele) then 'insertion') and (effect == 'frameshift_variant')) THEN 'frameshift-non-synonymous',
           WHEN ((length(ref) > length(allele) then 'insertion') and (effect == 'disruptive_inframe_insertion')) THEN 'non-synonymous',
           WHEN ((length(ref) < length(allele)) and (effect == 'frameshift_variant')) THEN 'frameshift-non-synonymous',
           WHEN ((length(ref) < length(allele)) and (effect == 'disruptive_inframe_deletion')) THEN 'non-synonymous',
           ELSE 'synonymous' END AS ppc_rating
        FROM wgs_ilmn.vars_coding v
        WHERE v.chrom = '{chrom}'
        AND v.blk_pos = {pos};
        '''.format(chrom, pos)
        snpeff.run_query(add_ppc)

snpeff.run_query("compute stats wgs_ilmn.global_vars")