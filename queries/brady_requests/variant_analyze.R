#!/usr/bin/env Rscript

###############################
## Read in command line args ##
###############################

##Grab user args from command line
args=(commandArgs(TRUE))

##print help if both args not provided 
if(length(args) < 3) {
  cat("\n **Enter arguments for gene, db and kav:** \n")
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      variant_analyze.R annotate variants then test for zygosity and MIE
      
      Arguments:
      --gene  - Genes of interest, surrounded by quotes, comma-sep, no spaces
      --db    - Database to search in, ex.p7_ptb.illumina_variant
      --kav   - Kaviar allele frequency cutoff  %
      --help             - print help text
      
      Example:
      ./variant_analyze.R --gene='BRCA1' --db='p7_ptb.illumina_variant' --kav=10 \n\n")
 
  q(save="no")
}
 
## Parse arguments and coerce to list
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1


##load libraries, install if needed
list.of.packages <- c("RODBC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.fhcrc.org")
 
#######################
## Connect to impala ##
#######################
#load library
library(RODBC)
#connect using the DSN name you created on your machine
conn <- odbcConnect("Impala DSN")
print("Connected to impala.")

#################################
## Create query from user args ##
#################################
genes = args$genes
database = args$db
kav_freq = args$kav

print(genes)
print(database)
print(kav_freq)

# #create df of genes of interest, mark if wildcard
# args.df = as.data.frame(matrix(unlist(strsplit(args, ","))))
# colnames(args.df) = "gene"
# args.df$type = NA
# args.df[(grepl("%$", args.df$gene)),]$type = "wildcard"
# args.df[!(grepl("%$", args.df$gene)),]$type = "gene"
# 
# #create wildcard gene lists
# if (dim(args.df[(args.df$type == "wildcard"),])[1] > 0 ){
#   wildcards = paste("ens.gene_name LIKE '", paste(as.character(args.df[(args.df$type == "wildcard"),]$gene), 
#                                                   collapse= "' OR ens.gene_name LIKE '"), "'", sep="")
#   wildcards
# } else
#   wildcards = NULL
# 
# #create non-wildcard gene list
# if (dim(args.df[(args.df$type == "gene"),])[1] > 0){
#   genes = paste("WHERE ens.gene_name IN ('", paste(as.character(args.df[(args.df$type == "gene"),]$gene), collapse="','")
#                 , "')", sep="")
# } else
#   genes = NULL
# 
# 
# #if wildcard and genes, add statement to gene list
# if ((length(wildcards) > 0) & length(genes) > 0 ){
#   gene_list = paste(genes, "OR", wildcards, sep=" ")
# } else if 
#  (length(wildcards) > 0 & length(genes) < 1){
#   gene_list = paste("WHERE", wildcards)
# } else if (length(wildcards) < 1 & length(genes) > 0){
#   gene_list = genes
# } else 
#   print("No genes of interest found.")
# 
# #######################
# ##  Run impala query ##
# #######################
# 
# print(sqlQuery(conn, query))
# 
# 
# ben_query = paste0(" 
#   WITH ens AS (
#   SELECT DISTINCT chromosome as chr, start, stop, gene_name
#   FROM public_hg19.ensembl_genes
#       WHERE (gene_name IN (",  gene_list, 
#       " AND chromosome NOT LIKE 'H%'
# )
#       SELECT p.sample_id, p.qual, p.filter, k.id as rsID, (k.alle_freq * 100) as kav_pct, k.alle_cnt as kav_count,
#       gene_name, p.chr, p.pos, p.ref, p.alt, p.gt, 
#       (CASE  
#       WHEN SUBSTRING(p.sample_id, -2) = '01'
#       THEN 'M'
#       WHEN SUBSTRING(p.sample_id, -2) = '02'
#       THEN 'F'
#       WHEN SUBSTRING(p.sample_id, -2) = '03'
#       THEN 'NB'
#       END) as member, CONCAT(gene_name, ':', p.chr, ':', CAST(p.pos AS STRING)) as variant_id, 
#       CONCAT(p.chr, ':', CAST(p.pos AS STRING), ':', p.alt) as alt_id 
#       FROM
#       (SELECT DISTINCT p.sample_id, p.qual, p.filter, ens.gene_name, p.chr, p.pos, p.ref, p.alt, p.gt
#       FROM ens," paste(
#       WHERE p.chr = ens.chr
#       AND (p.pos >= ens.start AND p.pos <= ens.stop)
#       AND p.gt IS NOT NULL
#       ) AS p
#       LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
#       ON p.chr = k.chromosome
#       AND p.pos = k.pos
#       AND p.ref = k.ref
#       AND p.alt = k.alt
#       WHERE (k.alle_freq < .10 OR k.alle_freq IS NULL)
#       ")
# nchar(ben_query)
# 
# query =  strwrap(ben_query, width=nchar(ben_query), simplify=TRUE)
# query
# 
# WITH ens AS (
#   SELECT DISTINCT chromosome as chr, start, stop, gene_name
#   FROM public_hg19.ensembl_genes
#   WHERE (gene_name IN ( 'HADH', 'HADHA', 'HADHB', 'ACAA1',
#                         'ACAA2', 'EHHADH', 'ECHS1')
#          OR gene_name LIKE 'HSD17B%')
#   AND chromosome NOT LIKE "H%"
# )
# SELECT p.sample_id, p.qual, p.filter, k.id as rsID, (k.alle_freq * 100) as kav_pct, k.alle_cnt as kav_count,
# gene_name, p.chr, p.pos, p.ref, p.alt, p.gt, 
# (CASE  
#  WHEN SUBSTRING(p.sample_id, -2) = '01'
#  THEN 'M'
#  WHEN SUBSTRING(p.sample_id, -2) = '02'
#  THEN 'F'
#  WHEN SUBSTRING(p.sample_id, -2) = '03'
#  THEN 'NB'
#  END) as member, CONCAT(gene_name, ":", p.chr, ":", CAST(p.pos AS STRING)) as variant_id, 
# CONCAT(p.chr, ":", CAST(p.pos AS STRING), ":", p.alt) as alt_id 
# FROM
# (SELECT DISTINCT p.sample_id, p.qual, p.filter, ens.gene_name, p.chr, p.pos, p.ref, p.alt, p.gt
#  FROM ens, p7_ptb.itmi_102_puzzle p
#  WHERE p.chr = ens.chr
#  AND (p.pos >= ens.start AND p.pos <= ens.stop)
#  AND p.gt IS NOT NULL
# ) AS p
# LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
# ON p.chr = k.chromosome
# AND p.pos = k.pos
# AND p.ref = k.ref
# AND p.alt = k.alt
# WHERE (k.alle_freq < .10 OR k.alle_freq IS NULL)
# 
