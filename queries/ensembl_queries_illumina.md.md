Sample queries for ensembl gene annotations
===========================================

Find all variants in ensembl gene regions
-----------------------------------------

The following query can be used to determine what gene regions a variant
lies in and to create a subset for downstream annotation for finding
amino acid changes.

Update or remove the following variables, as needed:  
- vcf.filter: filter out variants that did not pass filtering by the
variant caller  
- vcf.qual: filter out variants below a certain quality score

    SELECT vcf.sample_id as vcf_sample_id, vcf.chromosome as vcf_chrom, vcf.position as vcf_pos,vcf.ref as vcf_ref, vcf.alt as vcf_alt, vcf.id as vcf_rsID, ens.start as ensembl_start, ens.stop as ensembl_end, ens.feature, ens.gene_name as ensembl_gene_name, ens.gene_id as ensembl_geneid, ens.gene_biotype as ensembl_gene_biotype, ens.transcript_name as ensembl_tx_name,ens.transcript_id as ensembl_trans_id, ens.exon_id as ensembl_exonid, ens.strand as ensembl_strand
    FROM p7_ptb.illumina_variant as vcf, public_hg19.ensembl_genes as ens  
    WHERE vcf.filter = "PASS" 
    WHERE vcf.qual > 100
    AND vcf.chromosome = ens.chromosome  
    AND vcf.position BETWEEN ens.start AND ens.stop  

This query is equivalent to:

    SELECT vcf.sample_id as vcf_sample_id, vcf.chromosome as vcf_chrom, vcf.position as vcf_pos,vcf.ref as vcf_ref, vcf.alt as vcf_alt, vcf.id as vcf_rsID, ens.start as ensembl_start, ens.stop as ensembl_end, ens.feature, ens.gene_name as ensembl_gene_name, ens.gene_id as ensembl_geneid, ens.gene_biotype as ensembl_gene_biotype, ens.transcript_name as ensembl_tx_name,ens.transcript_id as ensembl_trans_id, ens.exon_id as ensembl_exonid, ens.strand as ensembl_strand
    FROM p7_ptb.illumina_variant as vcf
    JOIN
    public_hg19.ensembl_genes as ens 
    ON vcf.chromosome = ens.chromosome  
    WHERE vcf.filter = "PASS"
    AND vcf.qual > 100
    AND vcf.position BETWEEN ens.start AND ens.stop 

Find specific genes in the ensembl gene subset created above
------------------------------------------------------------

Update or remove the following variables, as needed:  
- vcf.filter: filter out variants that did not pass filtering by the
variant caller  
- vcf.qual: filter out variants below a certain quality score  
- ens.gene\_name: enter each gene of interest as a list, with each gene
in quotes, comma-separated

    SELECT vcf.sample_id as vcf_sample_id, vcf.chromosome as vcf_chrom, vcf.position as vcf_pos,vcf.ref as vcf_ref, vcf.alt as vcf_alt, vcf.id as vcf_rsID, ens.start as ensembl_start, ens.stop as ensembl_end, ens.feature, ens.gene_name as ensembl_gene_name, ens.gene_id as ensembl_geneid, ens.gene_biotype as ensembl_gene_biotype, ens.transcript_name as ensembl_tx_name,ens.transcript_id as ensembl_trans_id, ens.exon_id as ensembl_exonid, ens.strand as ensembl_strand
    FROM p7_ptb.illumina_variant as vcf, public_hg19.ensembl_genes as ens  
    WHERE vcf.filter = "PASS"  
    AND vcf.qual > 100
    AND vcf.chromosome = ens.chromosome  
    AND vcf.position BETWEEN ens.start AND ens.stop 
    AND ens.gene_name IN ("RMRPP1","PPIAP13","NDST2","RP11-574K11.8","RPL39P25")

Find a specific gene in the ensebml gene subset
-----------------------------------------------

Update or remove the following variables, as needed:  
- vcf.filter: filter out variants that did not pass filtering by the
variant caller  
- vcf.qual: filter out variants below a certain quality score  
- ens.gene\_name: enter a gene of interest in quotes

    SELECT vcf.sample_id as vcf_sample_id, vcf.chromosome as vcf_chrom, vcf.position as vcf_pos,vcf.ref as vcf_ref, vcf.alt as vcf_alt, vcf.id as vcf_rsID, ens.start as ensembl_start, ens.stop as ensembl_end, ens.feature, ens.gene_name as ensembl_gene_name, ens.gene_id as ensembl_geneid, ens.gene_biotype as ensembl_gene_biotype, ens.transcript_name as ensembl_tx_name,ens.transcript_id as ensembl_trans_id, ens.exon_id as ensembl_exonid, ens.strand as ensembl_strand
    FROM p7_ptb.illumina_variant as vcf, public_hg19.ensembl_genes as ens  
    WHERE vcf.filter = "PASS"  
    AND vcf.qual > 100
    AND vcf.chromosome = ens.chromosome  
    AND vcf.position BETWEEN ens.start AND ens.stop 
    AND ens.gene_name IN ("RMRPP1")

Find variants in located in coding regions, exons or start and stop codons
--------------------------------------------------------------------------

Use this query to find all variants located in a region annotated by
ensembl as coding, exon, start or stop

Update or remove the following variables, as needed:  
- vcf.filter: filter out variants that did not pass filtering by the
variant caller  
- vcf.qual: filter out variants below a certain quality score  
- ens.feature: enter either "CDS", "exon", "start\_codon", "stop\_codon"
surrounded by quotes

    SELECT vcf.sample_id as vcf_sample_id, vcf.chromosome as vcf_chrom, vcf.position as vcf_pos,vcf.ref as vcf_ref, vcf.alt as vcf_alt, vcf.id as vcf_rsID, ens.start as ensembl_start, ens.stop as ensembl_end, ens.feature, ens.gene_name as ensembl_gene_name, ens.gene_id as ensembl_geneid, ens.gene_biotype as ensembl_gene_biotype, ens.transcript_name as ensembl_tx_name,ens.transcript_id as ensembl_trans_id, ens.exon_id as ensembl_exonid, ens.strand as ensembl_strand
    FROM p7_ptb.illumina_variant as vcf, public_hg19.ensembl_genes as ens  
    WHERE vcf.filter = "PASS"  
    AND vcf.qual > 100
    AND vcf.chromosome = ens.chromosome  
    AND vcf.position BETWEEN ens.start AND ens.stop 
    AND ens.feature = "CDS"

Find variant located in eiher one feature or another
----------------------------------------------------

Use this query to find all variants located in either one region, such
as coding, or another, such as exon.

Update or remove the following variables, as needed:  
- vcf.filter: filter out variants that did not pass filtering by the
variant caller  
- vcf.qual: filter out variants below a certain quality score  
- ens.feature: enter the fetures you are looking for "CDS", "exon",
"start\_codon", "stop\_codon", as shown below

    SELECT vcf.sample_id as vcf_sample_id, vcf.chromosome as vcf_chrom, vcf.position as vcf_pos,vcf.ref as vcf_ref, vcf.alt as vcf_alt, vcf.id as vcf_rsID, ens.start as ensembl_start, ens.stop as ensembl_end, ens.feature, ens.gene_name as ensembl_gene_name, ens.gene_id as ensembl_geneid, ens.gene_biotype as ensembl_gene_biotype, ens.transcript_name as ensembl_tx_name,ens.transcript_id as ensembl_trans_id, ens.exon_id as ensembl_exonid, ens.strand as ensembl_strand
    FROM p7_ptb.illumina_variant as vcf, public_hg19.ensembl_genes as ens  
    WHERE vcf.filter = "PASS"  
    AND vcf.qual > 100
    AND vcf.chromosome = ens.chromosome  
    AND vcf.position BETWEEN ens.start AND ens.stop 
    AND (ens.feature = "CDS" or ens.feature = "exon")

Find variants that are not located in ensembl gene regions
----------------------------------------------------------

This query can be used downstream to locate intergenic variants by
finding the two flanking genes, and distances between the variants and
the flanking genes.

Update or remove the following variables, as needed:  
- vcf.filter: filter out variants that did not pass filtering by the
variant caller  
- vcf.qual: filter out variants below a certain quality score

    SELECT vcf.* 
    FROM
      p7dev.illumina_test as vcf
      LEFT JOIN p7dev.ensembl_test as ens
        ON vcf.chromosome = ens.chromosome
        AND vcf.position BETWEEN ens.start AND ens.stop
    WHERE
      vcf.filter = "PASS" 
      AND vcf.qual > 100
      AND ens.chromosome IS NULL
