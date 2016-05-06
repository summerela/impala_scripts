-- create kav table
create table p7_product.kav_distinct
(
  pos      INT,
  ref STRING,
  alt STRING,
  kav_freq FLOAT,
  kav_source STRING
)
  PARTITIONED BY (chrom string, pos_block int);

-- insert kaviar subset into table by chromosome using bash
#for x in $(seq 1 22) X Y; do echo "$x"; nohup impala-shell -q "\
insert into table p7_product.kav_distinct partition (chrom, pos_block)
  select distinct pos, ref, alt, allele_freq as kav_freq,
    CASE when sources is NULL then 'dbSNP'
      else sources
      END as kav_source, chrom, cast(pos/1000000 as INT) as pos_block
  from p7_ref_grch37.kaviar_isb
  where chrom = '$x';
# "; done

-- compute stats on new table
compute stats p7_product.kav_distinct;

-- CREATE CLIN_DISTINCT
create table p7_product.clin_distinct
  (pos int,
  ref string,
  alt string,
  clin_sig string,
  clin_dbn string
)
  partitioned by (chrom string, pos_block int);

-- insert distinct subset into table
insert into table p7_product.clin_distinct partition (chrom, pos_block)
  select distinct pos, ref, alt, clin_sig, clin_dbn,
    CASE
      when chrom = 'MT' then 'M'
        else chrom
    END as chrom,
    cast(pos/1000000 as INT) as pos_block
  from p7_ref_grch37.clinvar;

-- compute stats on new table
compute stats p7_product.clin_distinct;

-- CREATE DBSNP SUBSET
create table p7_product.dbsnp_distinct
    (pos int,
    ref string,
    alt string,
    rs_id string,
    dbsnp_buildid string
    )
  partitioned by (chrom string, pos_block int);

-- insert subset into table
insert into p7_product.dbsnp_distinct partition (chrom, pos_block)
select pos, ref, alt, rs_id, dbsnpbuildid as dbsnp_build
  CASE
      when chrom = 'MT' then 'M'
        else chrom
    END as chrom,
  cast(pos/1000000 as INT) as pos_block
  from p7_ref_grch37.dbsnp;

-- compute stats on new table
compute stats p7_product.dbsnp_distinct;

-- CREATE DBNSPF_DISTINCT 
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

-- insert dbnsfp subset into table
#for x in $(seq 1 22) M MT X Y; do nohup impala-shell -q "\
INSERT INTO p7_product.dbnsfp_distinct PARTITION (chrom, pos_block)
    select distinct pos, ref, alt, cadd_raw, dann_score,
      interpro_domain, chrom, cast(pos/1000000 as int) as pos_block
    from p7_ref_grch37.dbnsfp_variant
    where chrom = '$x';
# " ; done

--- compute stats on new table
compute stats p7_product.dbnsfp_distinct;

-- create ens_distinct
create table p7_product.ens_distinct
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

INSERT INTO p7_product.ens_distinct partition (chrom, pos_block)
    select distinct start, stop, strand, gene_name, gene_id,
      transcript_name, transcript_id, exon_id as exon_name, exon_number,
      CASE when chrom = 'MT' then 'M'
        else chrom
          END as chrom,  cast(start/1000000 as int) as pos_block
    from p7_ref_grch37.ensembl_genes;

--- compute stats on new table
compute stats p7_product.ens_distinct;

-- create HGMD table partitioned by chrom and position
create table p7_product.hgmd_distinct
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
insert into table p7_product.hgmd_distinct partition (chrom, pos_block)
select pos, id as hgmd_id, ref, alt, var_class as hgmd_class,
  mut_type as hgmd_mut, dna as hgmd_dna, prot as hgmd_prot, phen as hgmd_phen,
  chrom, cast(pos/1000000 as int) as pos_block
  from p7_ref_grch37.hgmd;

-- compute stats on new table
compute stats p7_product.hgmd_distinct;