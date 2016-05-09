# noinspection SqlNoDataSourceInspectionForFile

-- TODO: store all tables as parquet
-- TODO: check which tables need M chrom

-- GLOBAL VARIANTS TABLE PIPELINE


--- create blank table
create table anno_grch37.kav_vars
   (pos int,
   ref string,
   alt string,
   kav_freq float,
   kav_source string)
partitioned by (chrom string, pos_block int);

--- insert variants into anno_grch37.kav_vars and annotate with kaviar frequency and source

#for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into anno_grch37.kav_vars partition (chrom, pos_block)
SELECT CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt,
       t0.kav_freq,
       t0.kav_source,
       coalesce(t0.chrom, t1.chrom) AS chrom,
       CAST(coalesce(t0.pos, t1.pos)/1000000 AS int) as pos_block
FROM anno_grch37.kav_distinct t0
FULL OUTER JOIN anno_grch37.distinct_vars t1
    ON t0.chrom = t1.chrom
    AND t0.pos = t1.pos
    AND t0.ref = t1.ref
    AND t0.alt = t1.alt
where t0.chrom = '$x'
and t0.pos_block = $y;
# "; done

--- compute stats on new table
compute stats anno_grch37.kav_vars;

-- CREATE CLIN_KAV_VARS (illumina + kaviar + clinvar)
create table anno_grch37.clin_vars
   (pos int,
   ref string,
   alt string,
   kav_freq float,
   kav_source string,
   clin_sig string,
   clin_dbn string)
partitioned by (chrom string, pos_block int);

-- insert variants into table
#for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table anno_grch37.clin_vars partition (chrom, pos_block)
SELECT CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt,
        t0.kav_freq,
        t0.kav_source,
        t1.clin_sig,
        t1.clin_dbn,
        coalesce(t0.chrom, t1.chrom) AS chrom,
        coalesce(t0.pos_block, t1.pos_block) as pos_block
FROM anno_grch37.kav_vars t0
FULL OUTER JOIN anno_grch37.clin_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt
WHERE t0.chrom = '$x'
AND t0.pos_block = $y;
# "; done

-- compute stats on new table
compute stats anno_grch37.clin_vars;

-- 694405468 rows in table

-- CREATE DBSNP_VARS (illumini, cgi, kaviar, clinvar + dbSNP)
-- columns = chrom, pos, ref, alt, kav_freq, kav_source, clin_sig, clin_dbn, rs_id, dbsnpbuildid as dbnsp_build

-- add rsID and variants from dbsnp to create all_variants table
create table anno_grch37.dbsnp_vars
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
partitioned by (chrom string, pos_block int);

--- insert variants into partitioned table
#for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into anno_grch37.dbsnp_vars partition (chrom, pos_block)
SELECT CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt,
        t1.rs_id,
        t0.clin_sig,
        t0.clin_dbn,
        t0.kav_freq,
        t0.kav_source,
        t1.dbsnp_buildid as dbsnp_build,
        coalesce(t0.chrom, t1.chrom) AS chrom,
        coalesce(t0.pos_block, t1.pos_block) AS pos_block
FROM anno_grch37.clin_vars t0
FULL OUTER JOIN anno_grch37.dbsnp_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt
where t0.chrom = '$x'
and t0.pos_block = $y;
# " ; done

--- compute stats on new table
compute stats anno_grch37.dbsnp_vars;

-- 698,889,772 rows in table, all are distinct

-- CREATE ALL_VARS (all_vars + dbnsfp)
-- partitions = chrom, pos_block

create table anno_grch37.all_vars
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    transcript_name string,
    transcript_id string,
    exon_name string,
    exon_number int,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    cadd_raw float,
    dann_score float,
    interpro_domain string
    )
    partitioned by (chrom string, pos_block int);

#for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table anno_grch37.all_vars partition (chrom, pos_block)
SELECT DISTINCT COALESCE(ens.pos, d.pos) as pos, COALESCE(ens.ref, d.ref) as ref, COALESCE(ens.alt, d.alt) as alt,
ens.rs_id, ens.strand, ens.gene_name, ens.gene_id, ens.transcript_name, ens.transcript_id, ens.exon_name, ens.exon_number,
  ens.clin_sig, ens.clin_dbn, ens.kav_freq, ens.kav_source,
  ens.dbsnp_build, d.cadd_raw, d.dann_score, d.interpro_domain,
  COALESCE(ens.chrom,d.chrom) as chrom, ens.pos_block
from anno_grch37.ens_vars ens
full outer join anno_grch37.dbnsfp_distinct d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where ens.chrom = '$x'
and ens.pos_block = $y
and d.chrom = '$x'
and d.pos_block = $y;
# " ; done

-- compute stats on table
compute stats anno_grch37.all_vars;


-- ANNOTATE ALL_VARIANTS TABLE TO CREATE GLOBAL VARIANTS TABLE

-- create ens_distinct
create table anno_grch37.ens_distinct
(
  start INT,
  stop INT,
  strand STRING,
  gene_name STRING,
  gene_id STRING,
  transcript_name string,
  transcript_id string,
  exon_name string,
  exon_number int
)
  PARTITIONED BY (chrom string, pos_block int);

INSERT INTO anno_grch37.ens_distinct partition (chrom, pos_block)
    select distinct start, stop, strand, gene_name, gene_id,
      transcript_name, transcript_id, exon_id as exon_name, exon_number,
      CASE when chrom = 'MT' then 'M'
        else chrom
          END as chrom,  cast(start/1000000 as int) as pos_block
    from anno_grch37.ensembl_genes;

--- compute stats on new table
compute stats anno_grch37.ens_distinct;

-- 2243527 rows in table: PASS

-- CREATE ENS_VARS (all_vars + ens_distinct)
-- paritions = chrom, pos_block

create table anno_grch37.ens_vars
    (pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    transcript_name string,
    transcript_id string,
    exon_name string,
    exon_number int,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    cadd_raw float,
    dann_score float,
    interpro_domain string
    )
    partitioned by (chrom string, pos_block int);

-- inserted each chromosome into partitioned table as follows
#for x in $(seq 5 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table anno_grch37.ens_vars partition (chrom, pos_block)
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
  t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id, t1.exon_name, t1.exon_number,
  t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
  t0.cadd_raw, t0.dann_score, t0.interpro_domain, t0.chrom, t0.pos_block
FROM anno_grch37.all_vars t0
LEFT JOIN anno_grch37.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE t0.chrom = '$x'
AND t0.pos_block = $y
AND t1.chrom = '$x'
and t1.pos_block = $y;
# " ; done ; done

-- compute stats on ens_vars
compute stats anno_grch37.ens_vars;

-- check that row counts are as expected by comparing results with all_vars table
-- 776,862,245 rows in all_vars table
-- 775,742,611 distinct rows in table, will use distinct clause on next table merge


--- ANNOTATE VARIANTS WITH HGMD ANNOTATIONS

-- create HGMD table partitioned by chrom and position
create table anno_grch37.hgmd_dist
(
  pos int,
  hgmd_id string,
  ref string,
  alt string,
  hgmd_class string,
  hgmd_mut string,
  hgmd_dna string,
  hgmd_prot string,
  hgmd_phen string
)
partitioned by (chrom string, pos_block int);

-- insert subset of hgmd columns into table
insert into table anno_grch37.hgmd_dist partition (chrom, pos_block)
select pos, id as hgmd_id, ref, alt, var_class as hgmd_class,
  mut_type as hgmd_mut, dna as hgmd_dna, prot as hgmd_prot, phen as hgmd_phen,
  chrom, cast(pos/1000000 as int) as pos_block
  from anno_grch37.hgmd;

-- compute stats on new table
compute stats anno_grch37.hgmd_dist;

-- join hgmd with dbsnpf_vars to create hgmd_vars
create table anno_grch37.hgmd_vars
(
    pos int,
    ref string,
    alt string,
    rs_id string,
    strand string,
    gene_name string,
    gene_id string,
    transcript_name string,
    transcript_id string,
    exon_name string,
    exon_number int,
    clin_sig string,
    clin_dbn string,
    kav_freq float,
    kav_source string,
    dbsnp_build int,
    cadd_raw float,
    dann_score float,
    interpro_domain string,
    hgmd_id string,
    hgmd_class string,
    hgmd_mut string,
    hgmd_dna string,
    hgmd_prot string,
    hgmd_phen string
)
partitioned by (chrom string, pos_block int)


-- left join dbnsfp_vars with hgmd annotations
#for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table anno_grch37.hgmd_vars partition (chrom, pos_block)
select d.pos, d.ref, d.alt, d.rs_id, d.strand, d.gene_name,d.gene_id, d.transcript_name,
    d.transcript_id, d.exon_name, d.exon_number, d.clin_sig, d.clin_dbn, d.kav_freq,
    d.kav_source, d.dbsnp_build, d.cadd_raw, d.dann_score, d.interpro_domain,
    h.hgmd_id, h.hgmd_class, h.hgmd_mut, h.hgmd_dna, h.hgmd_prot, h.hgmd_phen, d.chrom, d.pos_block
from anno_grch37.dbnsfp_vars d
left join anno_grch37.hgmd_dist h
    on d.chrom = h.chrom
    and d.pos = h.pos
    and d.ref = h.ref
    and d.alt = h.alt
where d.chrom = '$x'
and d.pos_block = $y;
#"; done; done

--- ANNOTATE VARIANTS WITH CODING CONSEQUENCES

--- create coding consequence table partitioned by chrom and position
create table anno_grch37.coding_partitioned
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
insert into table anno_grch37.coding_partitioned partition (chrom, pos_block)
  select pos, id, ref, alt, gene as gene_name, gene_id, effect,
    impact, feature, feature_id, biotype, rank, hgvs_c, hgvs_p,
    chrom, cast(pos/1000000 as int) as pos_block
  from anno_grch37.all_coding
  where chrom = '$x';
#" ; done

compute stats anno_grch37.coding_partitioned;
--- compared total rows and rows in each chromosome to ensure all rows were retained

--- create the final table in the pipeline, global_vars containing all basic annotations
CREATE TABLE anno_grch37.global_vars
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

#for x in $X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table anno_grch37.global_vars partition (chrom, pos_block, clin_sig, kav_rank, impact)
  WITH t1 as (
      SELECT *
      from anno_grch37.dbnsfp_vars
      where chrom = '$x'
      and pos_block = $y),
  t2 as (
      SELECT *
      from anno_grch37.coding_partitioned
      where chrom = '$x'
      and pos_block = $y)

SELECT t1.pos, t1.ref, t1.alt, t1.rs_id, t1.strand, t1.gene_name, t1.gene_id, t1.clin_dbn, t1.kav_freq, t1.kav_source,
  t1.dbsnp_build, t1.cadd_raw, t1.dann_score, t1.interpro_domain, t2.effect, t2.feature, t2.feature_id,
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

--- duplications were induced in table due to having to repeat queries to ensure no variants were missed when server crashed
--- creating new table with distinct variants and transforming MT chromosome to M
CREATE TABLE anno_grch37.global_variants
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

#for x in $(seq 1 22) X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into table anno_grch37.global_variants partition (chrom, pos_block, clin_sig, kav_rank, impact)
select distinct * from anno_grch37.global_vars where chrom = '$x' and pos_block = $y
# " ; done ; done

# nohup impala-shell -q "\
insert into table anno_grch37.global_variants partition (chrom, pos_block, clin_sig, kav_rank, impact)
select pos, ref, alt, rs_id, strand, gene_name, gene_id, clin_dbn, kav_freq, kav_source, dbsnp_build, cadd_raw, dann_score, interpro_domain,
EFFECT, FEATURE, FEATUREID, BIOTYPE, RANK, HGVS_C, HGVS_P, 'M' as chrom, pos_block, clin_sig, kav_rank, impact
  from anno_grch37.global_vars
  where chrom = 'MT'
# "

--- compute stats on new table
compute stats global_variants;


--- run unit tests

-- clean up temp tables
drop table kav_distinct;
drop table clin_distinct;
drop table dbsnp_distinct;
drop table clin_kav_distinct;



-- TODO: add distinct RCF variants here when available
-- TODO: reassess after normalization





