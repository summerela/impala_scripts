-- GLOBAL VARIANTS TABLE PIPELINE

-- SUBSET TABLES TO BE JOINED FOR MORE EFFICIENT QUERIES

-- create distinct subset of kaviar table 
create table kav_distinct as (
  select distinct chrom, pos, stop, ref, alt, 
  allele_freq as kav_freq, sources as kav_source
  from p7_ref_grch37.kaviar
  );

-- create distinct subset of clinvar table

create table clin_distinct as (
  select distinct chromosome as chrom, pos, ref, alt, clnsig as clin_sig, clndbn as clin_dbn
  from p7_ref_hg19.clinvar
  );

-- create distinct subset of dbsnp table
  create table dbsnp_distinct as (
  select distinct chrom, pos, ref, alt, 
  rs_id, dbsnpbuildid, vc
  from p7_ref_grch37.dbsnp
  );

-- add distinct RCF variants here when available

-- JOIN ALL VARIANT TABLES TO CREATE ALL_VARIANTS TABLE
 
-- outer join kaviar and clinvar tables  
-- TODO: reassess after normalization
create table training.clin_kav_distinct as
SELECT DISTINCT coalesce(t0.chrom, t1.chrom) AS chrom,
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.ref) AS alt, t1.clin_sig,
       t1.clin_dbn, kav_freq, kav_source
FROM training.kav_distinct t0
FULL OUTER JOIN training.clin_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt;

-- add rsID and variants from dbsnp to create all_variants table
create table all_variants as 
SELECT coalesce(t0.chrom, t1.chrom) AS chrom,
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.ref) AS alt, t1.rs_id, t0.clin_sig,
       t0.clin_dbn, t0.kav_freq, t0.kav_source,
       t1.dbsnpbuildid AS dbsnp_build, t1.vc AS var_type
FROM training.clin_kav_distinct t0
FULL OUTER JOIN training.dbsnp_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt
WHERE t0.chrom = '1';

-- performed insert for each chromosome as follows
insert into table training.all_variants  
SELECT distinct coalesce(t0.chrom, t1.chrom) AS chrom,
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.ref) AS alt, t1.rs_id, t0.clin_sig,
       t0.clin_dbn, t0.kav_freq, t0.kav_source,
       t1.dbsnpbuildid AS dbsnp_build, t1.vc AS var_type
FROM training.clin_kav_distinct t0
FULL OUTER JOIN training.dbsnp_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt
WHERE t0.chrom = '2';

-- paritioned table for downstream analysis (do all in one step next time)
create table training.vars_partitioned
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

insert into training.vars_partitioned partition (chrom) 
select pos, ref, alt, rs_id, clin_sig, clin_dbn, kav_freq, kav_source,
dbsnp_build, var_type, chrom
from training.all_variants;

-- clean up temp tables
drop table kav_distinct;
drop table clin_distinct;
drop table dbsnp_distinct;
drop table clin_kav_distinct;

-- ANNOTATE ALL_VARIANTS TABLE TO CREATE GLOBAL VARIANTS TABLE

-- create distinct subset of ensembl table
create table ens_distinct as
select distinct chrom, start, stop, strand, gene_name, gene_id, transcript_name, transcript_id
from p7_ref_grch37.ensembl_genes;

-- create distinct subset of dbsnfp table
create table dbnsfp_distinct as (
  select distinct chrom, pos, ref, alt, cadd_raw, dann_score, interpro_domain
  from p7_ref_grch37.dbnsfp_variant
  );

-- add ensembl annotations to vars_partitioned to create ens_partitioned table
create table training.ens_partitioned
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
insert into table training.ens_partitioned partition (chrom)
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
    t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
    t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM training.vars_partitioned t0
LEFT JOIN training.ens_distinct t1
ON t0.chrom = t1.chrom
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '1';

-- add dbsnfp annotations to create final global_vars table
create table training.global_vars
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
    partitioned by (chrom string)

insert into table training.global_vars partition (chrom)
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from training.ens_partitioned ens
left join training.dbnsfp_distinct d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where ens.chrom = '1'
and d.chrom = '1';


-- clean up temp tables
drop table dbnsfp_distinct;
drop table ens_distinct;
drop table ens_partitioned;
drop table vars_partitioned;
-- drop table all_variants; drop this after using it create snpeff annotations