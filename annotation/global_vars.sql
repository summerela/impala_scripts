-- GLOBAL VARIANTS TABLE PIPELINE

-- SUBSET TABLES TO BE JOINED FOR MORE EFFICIENT QUERIES

-- create distinct subset of kaviar table 
create table p7_product.kav_distinct as (
  select distinct chrom, pos, stop, ref, alt, 
  allele_freq as kav_freq, sources as kav_source
  from p7_ref_grch37.kaviar_isb
  );

-- create distinct subset of clinvar table

create table p7_product.clin_distinct as (
  select distinct chrom, pos, ref, alt, clin_sig, clin_dbn
  from p7_ref_grch37.clinvar
  );

-- create distinct subset of dbsnp table
  create table p7_product.dbsnp_distinct as (
  select distinct chrom, pos, ref, alt, 
  rs_id, dbsnpbuildid, vc
  from p7_ref_grch37.dbsnp
  );

--- compute stats on new tables
compute stats p7_product.kav_distinct;
compute stats p7_product.clin_distinct;
compute stats p7_product.dbsnp_distinct;

-- add distinct RCF variants here when available

-- JOIN ALL VARIANT TABLES TO CREATE ALL_VARIANTS TABLE
-- TODO: reassess after normalization

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
insert into p7_product.clin_kav_distinct partition (chrom) 
SELECT DISTINCT 
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt, 
       t1.clin_sig,
       t1.clin_dbn, 
       t0.kav_freq, 
       t0.kav_source, 
       coalesce(t0.chrom, t1.chrom) AS chrom
FROM p7_product.kav_distinct t0
FULL OUTER JOIN p7_product.clin_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos = t1.pos 
    AND t0.ref = t1.ref 
    AND t0.alt = t1.alt
WHERE t0.chrom = '1'; 

--- repeat insert for each chromosome (only 2 shown here)
insert into p7_product.clin_kav_distinct partition (chrom) 
SELECT DISTINCT 
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt, 
       t1.clin_sig,
       t1.clin_dbn, 
       t0.kav_freq, 
       t0.kav_source, 
       coalesce(t0.chrom, t1.chrom) AS chrom
FROM p7_product.kav_distinct t0
FULL OUTER JOIN p7_product.clin_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos = t1.pos 
    AND t0.ref = t1.ref 
    AND t0.alt = t1.alt
WHERE t0.chrom = '2';

--- compute stats on new table
compute stats p7_product.clin_kav_distinct; 

-- add rsID and variants from dbsnp to create all_variants table
-- paritioned table for downstream analysis, did not do at creating due to coalesce of common columns
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

--- insert variants into partitioned table in pairs of three (fastest)
insert into table p7_product.all_vars partition(chrom)
SELECT DISTINCT 
       CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.alt) AS alt, 
        t1.rs_id, 
        t0.clin_sig,
        t0.clin_dbn, 
        t0.kav_freq, 
        t0.kav_source,
        t1.dbsnpbuildid AS dbsnp_build, 
        t1.vc AS var_type,
        coalesce(t0.chrom, t1.chrom) AS chrom
FROM p7_product.clin_kav_distinct t0
FULL OUTER JOIN p7_product.dbsnp_distinct t1
    ON t0.chrom = t1.chrom AND
       t0.pos = t1.pos AND
       t0.ref = t1.ref AND
       t0.alt = t1.alt
WHERE (t0.chrom = '1' or t0.chrom = '2' or t0.chrom = '3')
AND (t1.chrom = '1' or t1.chrom = '2' or t1.chrom = '3');

--- compute stats on new table
compute stats p7_product.all_vars; 

-- clean up temp tables
drop table kav_distinct;
drop table clin_distinct;
drop table dbsnp_distinct;
drop table clin_kav_distinct;

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
SELECT t0.pos, t0.ref, t0.alt, t0.rs_id, t1.strand,
    t1.gene_name, t1.gene_id, t1.transcript_name, t1.transcript_id,
    t0.clin_sig, t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build,
    t0.var_type, t0.chrom
FROM p7_product.all_vars t0
LEFT JOIN p7_product.ens_distinct t1
    ON t0.chrom = t1.chrom
    AND (t0.pos BETWEEN t1.start and t1.stop)
WHERE (t0.chrom = '1' or t0.chrom = '2' or t0.chrom = '3')
AND (t1.chrom = '1' or t1.chrom = '2' or t1.chrom = '3');

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
    partitioned by (chrom string)

insert into table p7_product.global_variants partition (chrom)
SELECT distinct ens.pos, ens.ref, ens.alt, ens.rs_id, ens.strand, ens.gene_name,
    ens.gene_id, ens.transcript_name, ens.transcript_id, ens.clin_sig,
    ens.clin_dbn, ens.kav_freq, ens.kav_source, ens.dbsnp_build, ens.var_type,
    d.cadd_raw, d.dann_score, d.interpro_domain, ens.chrom
from p7_product.ens_partitioned ens
left join p7_product.dbnsfp_distinct d
    on ens.chrom = d.chrom
    and ens.pos = d.pos
    and ens.ref = d.ref
    and ens.alt  = d.alt
where (ens.chrom = '1' or ens.chrom = '2' or ens.chrom = '3')
and (d.chrom = '1' or d.chrom = '2' or d.chrom = '3');


-- clean up temp tables
drop table dbnsfp_distinct;
drop table ens_distinct;
drop table ens_partitioned;
drop table vars_partitioned;

--- ANNOTATE VARIANTS WITH CODING CONSEQUENCES

--- download chrom, pos, rs_id as id, ref, alt from all_vars table
#nohup impala-shell -o all_vars -B -q 'select distinct chrom, pos, rs_id as id, ref, alt from p7_product.all_vars;' 

--- run snpeff
--- /tools/java/jdk1.7/bin/java -Xmx32g -jar /users/selasady/tools/snpEff/snpEff.jar -t -v -noStats GRCh37.74 all_vars > allvars_snpeff.vcf

-- turn snpeff output into a table
-- cat allvars_snpeff.vcf | /users/selasady/tools/snpEff/scripts/vcfEffOnePerLine.pl | /tools/java/jdk1.7/bin/java -Xmx32g -jar /users/selasady/tools/snpEff/SnpSift.jar extractFields \
--     - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].EFFECT" "ANN[*].IMPACT" \
--     "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" \
--     "ANN[*].HGVS_C" "ANN[*].HGVS_P" > dataset.tsv

--- make directory for storing file on hdfs
--- hdfs dfs -mkdir /user/selasady/global_vars/
--- hdfs dfs -put ./dataset.tsv /user/selasady/global_vars/

--- sql used to create table
CREATE TABLE p7_product.allvars_snpeff
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
COMMENT "snpeff output for all_vars table"
ROW FORMAT DELIMITED
    FIELDS TERMINATED BY '\t'

--- join snpeff coding consequences with global_variants table to create global_vars

--- compute stats on new table
compute stats dataset_snpeff;

-- drop table all_vars; 








