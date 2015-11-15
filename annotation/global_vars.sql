-- create global variants table 

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

  -- outer join kaviar and clinvar tables for SNV's 
  create table clin_kav_distinct as
  SELECT DISTINCT coalesce(t0.`chrom`, t1.`chrom`) AS `chrom`,
       CAST(coalesce(t0.`pos`, t1.`pos`) AS int) AS `pos`,
       coalesce(t0.`ref`, t1.`ref`) AS `ref`,
       coalesce(t0.`alt`, t1.`ref`) AS `alt`, t1.`clin_sig`,
       t1.`clin_dbn`, t0.`allele_freq` AS `kav_freq`,
       t0.`sources` AS `kav_source`
FROM training.`kav_distinct` t0
FULL OUTER JOIN training.`clin_distinct` t1
    ON t0.`chrom` = t1.`chrom` AND
       t0.`pos` = t1.`pos` AND
       t0.`ref` = t1.`ref` AND
       t0.`alt` = t1.`alt`
WHERE length(t0.`alt`) = 1

--  outer join kaviar and clinvar tables for indels and insert into clin_kav_distinct
insert into clin_kav_distinct 
SELECT DISTINCT coalesce(t0.`chrom`, t1.`chrom`) AS `chrom`,
       CAST(coalesce(t0.`pos`, t1.`pos`) AS int) AS `pos`,
       coalesce(t0.`ref`, t1.`ref`) AS `ref`,
       coalesce(t0.`alt`, t1.`ref`) AS `alt`, t1.`clin_sig`,
       t1.`clin_dbn`, t0.`allele_freq` AS `kav_freq`,
       t0.`sources` AS `kav_source`
FROM training.`kav_distinct` t0
FULL OUTER JOIN training.`clin_distinct` t1
    ON t0.`chrom` = t1.`chrom` AND
       t0.`pos` = t1.`pos` AND
       t0.`ref` = t1.`ref` AND
       t0.`alt` = t1.`alt`
WHERE length(t0.`alt`) != 1

-- add rsID and variants from dbsnp to create all_variants table
create table all_variants as 
SELECT coalesce(t0.`chrom`, t1.`chrom`) AS `chrom`,
       CAST(coalesce(t0.`pos`, t1.`pos`) AS int) AS `pos`,
       coalesce(t0.`ref`, t1.`ref`) AS `ref`,
       coalesce(t0.`alt`, t1.`ref`) AS `alt`, t1.`rs_id`, t0.`clin_sig`,
       t0.`clin_dbn`, t0.`kav_freq`, t0.`kav_source`,
       t1.`dbsnpbuildid` AS `dbsnp_build`, t1.`vc` AS `var_type`
FROM training.`clin_kav_distinct` t0
FULL OUTER JOIN training.`dbsnp_distinct` t1
    ON t0.`chrom` = t1.`chrom` AND
       t0.`pos` = t1.`pos` AND
       t0.`ref` = t1.`ref` AND
       t0.`alt` = t1.`alt`


-- create distinct subset of ensembl table
create table ens_distinct as
select distinct chrom, start, stop, strand, gene_name, gene_id
from p7_ref_grch37.ensembl_genes

-- create distinct subset of dbsnfp table
create table dbnsfp_distinct as (
  select distinct chrom, pos, ref, alt, cadd_raw, dann_score, interpro_domain
  from p7_ref_grch37.dbnsfp_variant
  )

  -- add ensembl annotations to all_variants to create vars_ens table
SET COMPRESSION_CODEC=bzip2;

create table vars_ens as 
SELECT t0.chrom, t0.pos, t0.ref,
       t0.alt, t0.rs_id, t0.clin_sig, t0.clin_dbn, t0.kav_freq,
       t0.kav_source, t0.dbsnp_build, t0.var_type,    
       t1.gene_name AS ens_gene, t1.gene_id AS ens_geneid,
       t1.strand
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos >= t1.start 
    AND t0.pos <= t1.stop

-- add dbsnfp annotations to create final global_vars table

SELECT t0.*, t1.cadd_raw, t1.dann_score, t1.interpro_domain
FROM training.vars_ens t0
  LEFT JOIN training.dbnsfp_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt

-- clean up temp tables
drop table kav_distinct;
drop table clin_distinct;
drop table dbnsfp_distinct;
drop table dbsnp_distinct;
drop table clin_kav_distinct;
drop table ens_distinct;
-- drop table all_variants; drop this after using it create snpeff annotations