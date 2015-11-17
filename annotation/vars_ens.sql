insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '2';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '3';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '4';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '5';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '6';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '7';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '8';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '9';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '10';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '11';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '12';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '13';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '14';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '15';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '16';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '17';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '18';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '19';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '20';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '21';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = '22';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = 'MT';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = 'X';

insert into table training.vars 
SELECT t0.chrom, t0.pos, t0.ref, t0.alt,
    t0.rs_id, t1.strand, t1.gene_name, t1.gene_id,
    t1.transcript_name, t1.transcript_id, t0.clin_sig,
    t0.clin_dbn, t0.kav_freq, t0.kav_source, t0.dbsnp_build, t0.var_type
FROM training.all_variants t0
LEFT JOIN training.ens_distinct t1
    ON t0.chrom = t1.chrom 
    AND t0.pos BETWEEN t1.start and t1.stop
WHERE t0.chrom = 'Y';
