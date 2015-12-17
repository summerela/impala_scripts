# noinspection SqlNoDataSourceInspectionForFile
-- GLOBAL VARIANTS TABLE PIPELINE

-- SUBSET TABLES TO BE JOINED FOR MORE EFFICIENT QUERIES

-- create distinct subset of kaviar table

create table p7_product.kav_distinct
(
  pos      INT,
  stop     INT,
  ref STRING,
  alt STRING,
  kav_freq FLOAT,
  kav_source STRING
)
  PARTITIONED BY (chrom string);

insert into table p7_product.kav_distinct partition (chrom)
  select distinct pos, stop, ref, alt,
  allele_freq as kav_freq, sources as kav_source, chrom
  from p7_ref_grch37.kaviar_isb;

compute stats p7_product.kav_distinct;

-- create distinct subset of clinvar table

create table p7_product.clin_distinct
  (pos int,
  ref string,
  alt string,
  clin_sig string,
  clin_dbn string
)
  partitioned by (chrom string);

insert into table p7_product.clin_distinct partition (chrom)
  select distinct pos, ref, alt, clin_sig, clin_dbn, chrom
  from p7_ref_grch37.clinvar;

compute stats p7_product.clin_distinct;

-- create distinct subset of dbsnp table

create table p7_product.dbsnp_distinct
    (pos int,
    ref string,
    alt string,
    rs_id string,
    dbsnp_buildid string,
    var_type string
    )
  partitioned by (chrom string)

insert into p7_product.dbsnp_distinct partition (chrom)
select distinct pos, ref, alt,
  rs_id, dbsnpbuildid as dbsnp_buildid, vc as var_type, chrom
  from p7_ref_grch37.dbsnp;

compute stats p7_product.dbsnp_distinct;

# select count used to compare rows in all tables to original

-- TODO: add distinct RCF variants here when available
-- TODO: reassess after normalization

-- JOIN ALL VARIANT TABLES TO CREATE ALL_VARIANTS TABLE


--- create blank partitioned table
create table p7_product.clin_kav_distinct
   (pos int,
   ref string,
   alt string,
   clin_sig string,
   clin_dbn string,
   kav_freq float,
   kav_source string)
partitioned by (chrom string) 

--- outer join kaviar and clinvar tables and insert into p7_product.clin_kav_distinct by chrom
# uncomment bash lines and run in shell
#for x in $(seq 1 22) M MT X Y; do echo "$x"; impala-shell -q "\
insert into p7_product.clin_kav_distinct partition (chrom)
WITH t0 as
(select * from p7_product.kav_distinct where chrom = '$x'),
t1 as
(select * from p7_product.clin_distinct where chrom = '$x')
SELECT CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt,
       t1.clin_sig,
       t1.clin_dbn,
       t0.kav_freq,
       t0.kav_source,
       coalesce(t0.chrom, t1.chrom) AS chrom
FROM t0
FULL OUTER JOIN t1
    ON t0.chrom = t1.chrom
    AND t0.pos = t1.pos
    AND t0.ref = t1.ref
    AND t0.alt = t1.alt;
# "; done

--- compute stats on new table
compute stats p7_product.clin_kav_distinct; 

-- add rsID and variants from dbsnp to create all_variants table
create table p7_product.all_vars
 (pos int,
    ref string,
    alt string,
    rs_id string,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build string,
    var_type string)
partitioned by (chrom string);

--- insert variants into partitioned table
#for x in $(seq 1 22) M MT X Y; do echo "$x"; nohup impala-shell -q "\
insert into p7_product.all_vars partition (chrom)
WITH t0 as
(select * from p7_product.clin_kav_distinct where chrom = '$x'),
t1 as
(select * from p7_product.dbsnp_distinct where chrom = '$x')
SELECT DISTINCT
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt,
        t1.rs_id,
        t0.clin_sig,
        t0.clin_dbn,
        t0.kav_freq,
        t0.kav_source,
        t1.dbsnp_buildid,
        t1.var_type,
        coalesce(t0.chrom, t1.chrom) AS chrom
FROM t0
FULL OUTER JOIN t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt;
# " ; done

--- compute stats on new table
compute stats p7_product.all_vars; 


-- ANNOTATE ALL_VARIANTS TABLE TO CREATE GLOBAL VARIANTS TABLE

-- create distinct subset of ensembl table
create table p7_product.ens_distinct as
select distinct chrom, start, stop, strand, gene_name, gene_id, transcript_name, transcript_id
from p7_ref_grch37.ensembl_genes;

--- compute stats on new table
compute stats p7_product.ens_distinct;

-- create distinct subset of dbsnfp table
create table p7_product.dbnsfp_distinct as (
  select distinct chrom, pos, ref, alt, cadd_raw, dann_score, interpro_domain
  from p7_ref_grch37.dbnsfp_variant
  );

--- compute stats on new table
compute stats p7_product.dbnsfp_distinct;

-- add ensembl annotations to vars_partitioned to create ens_partitioned table
create table p7_product.ens_partitioned
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    transcript_name string,
    transcript_id string,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build string,
    var_type string
    )
    partitioned by (chrom string)

-- inserted each chromosome into partitioned table as follows
insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '1'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '1';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '2'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '2';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '3'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '3';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '4'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '4';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '5'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '5';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '6'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '6';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '7'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '7';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '8'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '8';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '9'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '9';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '10'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '10';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '11'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '11';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '12'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '12';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '13'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '13';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '14'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '14';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '15'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '15';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '16'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '16';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '17'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '17';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '18'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '18';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '19'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '19';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '20'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '20';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '21'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '21';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = '22'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = '22';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = 'X'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = 'X';

insert into table p7_product.ens_partitioned partition (chrom)
with t0 as (
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id,  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
    FROM p7_product.all_vars t0
    WHERE t0.chrom = 'Y'
  )

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
   t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t1.chrom = 'Y';

--- partition dbsnfp_distinct table
create table dbsnfp_partitioned
  (pos int,
  ref string,
  alt string,
  cadd_raw float,
  dann_score float,
  interpro_domain string
  )
  PARTITIONED BY (chrom string);

--- add dbsnfp_distinct into dbnsfp_partitioned
INSERT INTO p7_product.dbnsfp_partitioned partition (chrom)
SELECT DISTINCT pos, ref, alt, cadd_raw, dann_score, interpro_domain, chrom
FROM p7_product.dbnsfp_distinct;

compute stats dbnsfp_partitioned;

-- add dbsnfp annotations to create final global_vars table
create table p7_product.global_variants
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    transcript_name string,
    transcript_id string,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build string,
    var_type string,
    cadd_raw float,
    dann_score float,
    interpro_domain string
    )
    partitioned by (chrom string);


insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '1'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '1';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '2'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '2';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '3'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '3';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '4'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '4';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '5'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '5';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '6'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '6';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '7'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '7';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '8'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '8';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '9'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '9';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '10'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '10';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '11'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '11';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '12'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '12';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '13'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '13';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '14'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '14';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '15'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '15';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '16'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '16';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '17'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '17';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '18'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '18';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '19'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '19';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '20'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '20';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '21'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '21';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = '22'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = '22';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = 'X'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = 'X';

insert into table p7_product.global_variants partition (chrom)
WITH ens as (
    SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    ens.chrom
    from p7_product.ens_partitioned ens
    where ens.chrom = 'Y'
  )
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from ens
left join p7_product.dbnsfp_partitioned d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where d.chrom = 'Y';

compute stats p7_product.global_variants;

--- unit test


-- clean up temp tables
drop table dbnsfp_distinct;
drop table ens_distinct;
drop table ens_partitioned;
drop table vars_partitioned;

--- ANNOTATE VARIANTS WITH CODING CONSEQUENCES

--- sql used to create table
CREATE TABLE p7_product.global_vars
(
  CHROM string ,
  POS int ,
  REF string ,
  ALT string ,
  GENE string ,
  GENEID string ,
  EFFECT string ,
  IMPACT string ,
  FEATURE string ,
  FEATUREID string ,
  BIOTYPE string ,
  RANK int ,
  HGVS_C string ,
  HGVS_P string ) 
COMMENT "snpeff output for global_vars table"
ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'

--- join snpeff coding consequences with global_variants table to create global_vars

--- compute stats on new table
compute stats dataset_snpeff;

-- clean up temp tables
drop table kav_distinct;
drop table clin_distinct;
drop table dbsnp_distinct;
drop table clin_kav_distinct;


-- drop table all_vars; 








