#!/usr/bin/env bash

## GLOBAL VARIANTS TABLE PIPELINE for Illumina variants ##

#########################
###    Assumptions    ###
#########################
# cgi_snpeff.py in same directory as gv_cgi.sh
# run on cluster with access to impala shell
# all variant tables have a var_id column
# all variant tables are left-normalized, 1 based coords

## see assumptions under cgi_snpeff.py


############################################################
# create table of all variants found by joining distinct ###
###   cgi variants with distinct reference variants      ###
############################################################
# create ilmn_vars empty table

impala-shell -q "create table wgs_cg.cg_vars
(
  var_id string,
  pos int,
  ref string,
  allele string)
partitioned by (chrom string, blk_pos int)
STORED AS PARQUET;"

# insert variants into wgs_cg.ilmn_vars
for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.cg_vars partition (chrom, blk_pos)
SELECT var_id, pos, ref, allele, chrom, blk_pos FROM wgs_cg.vcf_distinct WHERE chrom = '$x' AND blk_pos = $y
UNION
SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.dbnsfp_distinct WHERE chrom = '$x' AND blk_pos = $y
UNION
SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.kaviar_distinct WHERE chrom = '$x' AND blk_pos = $y
UNION
SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.clinvar_distinct WHERE chrom = '$x' AND blk_pos = $y
UNION
SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.dbsnp WHERE chrom = '$x' AND blk_pos = $y
UNION
SELECT var_id, pos, ref, alt as allele, chrom, blk_pos FROM anno_grch37.hgmd WHERE chrom = '$x' AND blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.cg_vars;"

###################
### RUN SNPEFF  ###
###################
# edit this script before running
python ./cg_snpeff.py &

##################################
### ADD REFERENCE ANNOTATIONS  ###
##################################

# rsID
impala-shell -q "
create table wgs_cg.vars_dbsnp
(
  var_id string,
  pos int,
  ref string,
  allele string,
  rs_id string,
  dbsnp_buildid string
)
partitioned by (chrom string, blk_pos int)
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_dbsnp partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, d.rs_id,
    d.dbsnpbuildid as dbsnp_buildid,
    v.chrom, v.blk_pos
FROM wgs_cg.cg_vars v
LEFT JOIN anno_grch37.dbsnp d
  ON v.var_id = d.var_id
WHERE v.chrom = '$x' and d.chrom = '$x'
AND v.blk_pos = $y and d.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_dbsnp;"

# kaviar
impala-shell -q "
create table wgs_cg.vars_kaviar
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_kaviar partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
    v.dbsnp_buildid, k.kav_freq, k.kav_source,
    v.chrom, v.blk_pos
FROM wgs_cg.vars_dbsnp v
LEFT JOIN anno_grch37.kaviar_distinct k
  ON v.var_id = k.var_id
WHERE v.chrom = '$x' and k.chrom = '$x'
AND v.blk_pos = $y and k.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_kaviar;"

# drop temp table from previous step
impala-shell -q "drop table wgs_cg.vars_dbsnp;:"

# clinvar
impala-shell -q "
create table wgs_cg.vars_clinvar
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_clinvar partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
    v.dbsnp_buildid, v.kav_freq, v.kav_source,
    c.clin_sig, c.clin_dbn, v.chrom, v.blk_pos
FROM wgs_cg.vars_kaviar v
LEFT JOIN anno_grch37.clinvar_distinct c
  ON v.var_id = c.var_id
WHERE v.chrom = '$x' and c.chrom = '$x'
AND v.blk_pos = $y and c.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_clinvar;"

# drop temp table from previous step
impala-shell -q "drop table wgs_cg.vars_kaviar;"

# hgmd
impala-shell -q "
create table wgs_cg.vars_hgmd
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_hgmd partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
    v.dbsnp_buildid, v.kav_freq, v.kav_source,
    v.clin_sig, v.clin_dbn, h.id as hgmd_id, h.var_class as hgmd_varclass,
    v.chrom, v.blk_pos
FROM wgs_cg.vars_clinvar v
LEFT JOIN anno_grch37.hgmd h
    ON v.var_id = h.var_id
WHERE v.chrom = '$x' and h.chrom = '$x'
AND v.blk_pos = $y and h.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_hgmd;"

# drop temp table from previous step
impala-shell -q "drop table wgs_cg.vars_clinvar;"

# dbnsfp cadd, dann, interpro
impala-shell -q "
create table wgs_cg.vars_dbnsfp
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_dbnsfp partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
    v.dbsnp_buildid, v.kav_freq, v.kav_source,
    v.clin_sig, v.clin_dbn, v.hgmd_id, v.hgmd_varclass,
    d.cadd_raw, d.dann_score, d.interpro_domain, v.chrom, v.blk_pos
FROM wgs_cg.vars_hgmd v
LEFT JOIN anno_grch37.dbnsfp_distinct d  ON v.var_id = d.var_id
WHERE v.chrom = '$x' and d.chrom = '$x'
AND v.blk_pos = $y and d.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_dbsnfp;"

# drop temp table from previous step
impala-shell -q "drop table wgs_cg.vars_hgmd;"

# ensembl gene, tx, exon
impala-shell -q "
create table wgs_cg.vars_ensembl
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_ensembl partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
    v.dbsnp_buildid, v.kav_freq, v.kav_source,
    v.clin_sig, v.clin_dbn, v.hgmd_id, v.hgmd_varclass,
    v.cadd_raw, v.dann_score, v.interpro_domain,
    e.strand, e.gene_name, e.gene_id, e.transcript_name, e.transcript_id,
    e.exon_name, e.exon_number, v.chrom, v.blk_pos
FROM wgs_cg.vars_dbsnfp v
LEFT JOIN anno_grch37.ensembl_distinct e
    ON v.chrom = e.chrom
    AND (v.pos BETWEEN e.pos and e.stop)
WHERE v.chrom = '$x' and e.chrom = '$x'
AND v.blk_pos = $y and e.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_ensembl;"

# drop temp table from previous step
impala-shell -q "drop table wgs_cg.vars_dbsnfp;"

################################
### ADD CODING CONSEQUNECES  ###
################################

# coding consequences
impala-shell -q "
create table wgs_cg.vars_coding
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.vars_coding partition (chrom, blk_pos)
SELECT v.var_id, v.pos, v.ref, v.allele, v.rs_id,
    v.dbsnp_buildid, v.kav_freq, v.kav_source,
    v.clin_sig, v.clin_dbn, v.hgmd_id, v.hgmd_varclass,
    v.cadd_raw, v.dann_score, v.interpro_domain,
    v.strand, v.gene_name, v.gene_id, v.transcript_name, v.transcript_id,
    v.exon_name, v.exon_number, c.effect, c.impact, c.feature, c.feature_id, c.biotype,
    c.rank, c.hgvs_c, c.hgvs_p, v.chrom, v.blk_pos
FROM wgs_cg.vars_ensembl v
LEFT JOIN anno_grch37.coding c
    ON v.chrom = c.chrom
    AND v.pos = c.pos
    AND v.ref = c.ref
    AND v.allele = c.allele
WHERE v.chrom = '$x' and c.chrom = '$x'
AND v.blk_pos = $y and c.blk_pos = $y;"; done; done

# compute stats on new table
impala-shell -q "compute stats wgs_cg.vars_coding;"

# drop temp table from previous step
impala-shell -q "drop table wgs_cg.vars_ensembl;"

###################################
### ANNOTATE WITH PPC NOTATION  ###
###################################

impala-shell -q "
create table wgs_cg.global_vars
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
 STORED AS PARQUET;"

for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "
insert into wgs_cg.global_vars partition (chrom, blk_pos)
SELECT v.*,
    CASE WHEN ((length(ref) == length(allele)) and (effect IN ('codings_sequence_variant', 'chromosome',
    'missense_variant', 'initator_codon_variant', 'stop_retained_variant', 'transcript_variant')) THEN 'non-synonymous',
    WHEN ((length(ref) > length(allele) then 'insertion') and (effect == 'frameshift_variant')) THEN 'frameshift-non-synonymous',
    WHEN ((length(ref) > length(allele) then 'insertion') and (effect == 'disruptive_inframe_insertion')) THEN 'non-synonymous',
    WHEN ((length(ref) < length(allele)) and (effect == 'frameshift_variant')) THEN 'frameshift-non-synonymous',
    WHEN ((length(ref) < length(allele)) and (effect == 'disruptive_inframe_deletion')) THEN 'non-synonymous',
    ELSE 'synonymous' END AS ppc_rating
FROM wgs_cg.vars_coding v
WHERE v.chrom = '$x'
AND v.blk_pos = $y;"; done; done
