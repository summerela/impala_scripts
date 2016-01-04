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

--- compared number of rows in each chromosome with number of rows in kav_distinct to ensure rows were retained

--- add additional positional partition to all_vars table to avoid exceeding memory on impala
create table p7_product.all_variants
( pos int,
 ref string,
 alt string,
 rs_id string,
 clin_sig string,
 clin_dbn string,
 kav_freq float,
 kav_source string,
 dbsnp_buildid int,
 var_type string
 )
 partitioned by (chrom string, pos_block int);

insert into p7_product.all_variants partition (chrom, pos_block)
select pos, ref, alt, rs_id, clin_sig, clin_dbn, kav_freq,
    kav_source, cast(dbsnp_build as int), var_type, chrom,
    cast(pos/1000000 as int) as pos_block
from p7_product.all_vars;

compute stats all_variants;

-- ANNOTATE ALL_VARIANTS TABLE TO CREATE GLOBAL VARIANTS TABLE

-- create distinct subset of ensembl table
create table p7_product.ens_distinct
(
  start      INT,
  stop     INT,
  strand STRING,
  gene_name STRING,
  gene_id STRING
)
  PARTITIONED BY (chrom string, pos_block int);

INSERT INTO p7_product.ens_distinct partition (chrom)
    select distinct start, stop, strand, gene_name, gene_id,
      chrom,  cast(start/1000000 as int) as pos_block
    from p7_ref_grch37.ensembl_genes;

--- compute stats on new table
compute stats p7_product.ens_distinct;

-- add ensembl annotations to vars_partitioned to create ens_partitioned table
create table p7_product.ens_partitioned
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    var_type string
    )
    partitioned by (chrom string, pos_block int);

-- inserted each chromosome into partitioned table as follows
#for x in $(seq 1 22) M MT X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table p7_product.ens_partitioned partition (chrom, pos_block)
with t0 as (
    SELECT * FROM p7_product.all_variants
    WHERE chrom = '$x'
    AND pos_block = $y),
t1 as (
  SELECT * FROM p7_product.ens_distinct
  WHERE chrom = '$x'
  AND pos_block = $y)

SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
   t1.gene_name, t1.gene_id, t0.clin_sig, t0.clin_dbn, t0.kav_freq,
  t0.kav_source, t0.dbsnp_buildid as buildid, t0.var_type, t0.chrom, cast(pos/1000000 as int) as pos_block
FROM t0
LEFT JOIN t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop);
# " ; done ; done

compute stats ens_partitioned;

--- check that row counts are as expected by comparing results with all_vars table
# for x in $(seq 1 22) M MT X Y; do impala-shell -q  \
# "select count(*) from p7_product.ens_partitioned where chrom = '$x';"; done

##for x in $(seq 1 22) M MT X Y; do impala-shell -q \
##"select count(*) from p7_product.all_vars where chrom = '$x';"; done

--- due to server issues, had to redo chromosomes 18 through MT, taking distcint of table
create table p7_product.ens_vars
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    var_type string
    )
    partitioned by (chrom string, pos_block int);

#for x in $(seq 1 22) M MT X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table p7_product.ens_vars partition (chrom, pos_block)
select distinct *
from p7_product.ens_partitioned
where chrom = '$x'
and pos_block = $y;
#" ; done; done


-- create distinct subset of dbsnfp table
create table p7_product.dbnsfp_distinct
(
  pos int,
  ref string,
  alt string,
  cadd_raw float,
  dann_score FLOAT,
  interpro_domain string
)
PARTITIONED BY (chrom string, pos_block int);

#for x in $(seq 1 22) M MT X Y; do nohup impala-shell -q "\
INSERT INTO p7_product.dbnsfp_distinct PARTITION (chrom, pos_block)
    select distinct pos, ref, alt, cadd_raw, dann_score,
      interpro_domain, chrom, cast(pos/1000000 as int) as pos_block
    from p7_ref_grch37.dbnsfp_variant
    where chrom = '$x';
# " ; done

--- compute stats on new table
compute stats p7_product.dbnsfp_distinct;

-- add dbsnfp annotations to create global_variants table
create table p7_product.dbnsfp_vars
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    var_type string,
    cadd_raw float,
    dann_score float,
    interpro_domain string
    )
    partitioned by (chrom string, pos_block int);

#for x in $(seq 1 22) M MT X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table p7_product.dbnsfp_vars partition (chrom, pos_block)
WITH ens as (
    SELECT *
    from p7_product.ens_vars
    where chrom = '$x'
    and pos_block = $y),
d as (
    SELECT *
    from p7_product.dbnsfp_distinct
    where chrom = '$x'
    and pos_block = $y)

SELECT ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
  ens.gene_id, ens.clin_sig, ens.clin_dbn, ens.kav_freq, ens.kav_source,
  ens.dbsnp_build, ens.var_type, d.cadd_raw, d.dann_score, d.interpro_domain,
  ens.chrom, ens.pos_block
from ens
left join d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt;
# " ; done

compute stats p7_product.dbnsfp_vars;

--- ANNOTATE VARIANTS WITH CODING CONSEQUENCES

--- create coding consequence table partitioned by chrom and position
create table p7_product.coding_partitioned
(
  pos int,
  id string,
  ref string,
  alt string,
  gene_name string,
  gene_id string,
  effect string,
  impact string,
  feature string,
  feature_id string,
  biotype string,
  rank int,
  hgvs_c string,
  hgvs_p string
)
partitioned by (chrom string, pos_block int);

#for x in $(seq 1 22) M MT X Y; do nohup impala-shell -q "\
insert into table p7_product.coding_partitioned partition (chrom, pos_block)
  select pos, id, ref, alt, gene as gene_name, gene_id, effect,
    impact, feature, feature_id, biotype, rank, hgvs_c, hgvs_p,
    chrom, cast(pos/1000000 as int) as pos_block
  from p7_product.all_coding
  where chrom = '$x';
#" ; done

compute stats p7_product.coding_partitioned;
--- compared total rows and rows in each chromosome to ensure all rows were retained

--- create the final table in the pipeline, global_vars containing all basic annotations
CREATE TABLE p7_product.global_vars
   (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    var_type string,
    cadd_raw float,
    dann_score float,
    interpro_domain string,
    EFFECT string,
    FEATURE string,
    FEATUREID string,
    BIOTYPE string,
    RANK int,
    HGVS_C string ,
    HGVS_P string )
partitioned by (chrom string, pos_block int, clin_sig string, kav_rank string, impact string)
COMMENT "Annotated table of all variants from Kaviar_ISB, ClinVar and dbSNP."

#for x in $(seq 16 22) M MT X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table p7_product.global_vars partition (chrom, pos_block, clin_sig, kav_rank, impact)
  WITH t1 as (
      SELECT *
      from p7_product.dbnsfp_vars
      where chrom = '$x'
      and pos_block = $y),
  t2 as (
      SELECT *
      from p7_product.coding_partitioned
      where chrom = '$x'
      and pos_block = $y)

SELECT t1.pos, t1.ref, t1.alt, t1.rs_id, t1.strand, t1.gene_name, t1.gene_id, t1.clin_dbn, t1.kav_freq, t1.kav_source,
  t1.dbsnp_build, t1.var_type, t1.cadd_raw, t1.dann_score, t1.interpro_domain, t2.effect, t2.feature, t2.feature_id,
  t2.biotype, t2.rank, t2.hgvs_c, t2.hgvs_p, t1.chrom, t1.pos_block, t1.clin_sig,
    (CASE
       when (t1.kav_freq < .03 ) then 'under3'
       when (t1.kav_freq > .03 AND t1.kav_freq <.05) then '3to5'
       else 'over5'
       END) AS kav_rank, t2.impact
FROM t1
LEFT JOIN t2
  ON t1.chrom = t2.chrom
  AND t1.pos = t2.pos
  AND t1.ref = t2.ref
  AND t1.alt = t2.alt;
# " ; done ; done

--- compute stats on new table
compute stats global_vars;

-- clean up temp tables
drop table kav_distinct;
drop table clin_distinct;
drop table dbsnp_distinct;
drop table clin_kav_distinct;









