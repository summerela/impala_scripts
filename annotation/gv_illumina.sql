# noinspection SqlNoDataSourceInspectionForFile

-- GLOBAL VARIANTS TABLE PIPELINE for Illumina variants
-- mix of bash and sql

-- join distinct illumina variants with Kaviar

--- create kav_vars empty table
create table anno_grch37.kav_vars
   (pos int,
   ref string,
   allele string,
   kav_freq float,
   kav_source string)
partitioned by (chrom string, pos_block int);

--- insert variants into anno_grch37.kav_vars and annotate with kaviar frequency and source

#for x in $(seq 1 22) M X Y; do for y in $(seq 0 249); do nohup impala-shell -q "\
insert into anno_grch37.kav_vars partition (chrom, pos_block)
-- coalesce required in case variants are found on only one side of the join, merges common columns
SELECT CAST(coalesce(t0.pos, t1.pos) AS int) AS pos,
       coalesce(t0.ref, t1.ref) AS ref,
       coalesce(t0.alt, t1.allele) AS allele,
       t0.kav_freq,
       t0.kav_source,
       coalesce(t0.chrom, t1.chrom) AS chrom,
       CAST(coalesce(t0.pos, t1.pos)/1000000 AS int) as pos_block
FROM anno_grch37.kav_distinct t0
FULL OUTER JOIN wgs_ilmn.test_vars t1
    ON t0.chrom = t1.chrom
    AND t0.pos = t1.pos
    AND t0.ref = t1.ref
    AND t0.alt = t1.allele
where t0.chrom = '$x'
and t0.pos_block = $y;
# "; done; done

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
       coalesce(t0.allele, t1.alt) AS allele,
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

