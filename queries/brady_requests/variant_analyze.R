#!/usr/bin/env Rscript

###############################
## Read in command line args ##
###############################

##Grab user args from command line
args=(commandArgs(TRUE))

##print help if both args not provided 
if(length(args) < 4) {
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
      --plat  - Specify if platform is illumina or cgi
      --help             - print help text
      
      Example:
      ./variant_analyze.R --gene='BRCA1' --db='p7_ptb.illumina_variant' --kav=10 --plat='illumina' \n\n")
  
  q(save="no")
}

## Parse arguments and coerce to list
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

#assign args to variable
gene_arg = argsL$gene
db_arg = argsL$db
kav_arg = argsL$kav
plat_arg = argsL$plat

##install library if needed
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
conn
#################################
## Create query from user args ##
#################################

gene_arg = "HADH,HADHA,HADHB,ACAA1,ACAA2,EHHADH,ECHS1,HSD17B%"
db_arg = "p7_ptb.itmi_102_puzzle"
kav_arg = "10"
plat_arg = "illumina"

#create df of genes of interest, mark if wildcard
genes.df = as.data.frame(matrix(unlist(strsplit(gene_arg, ","))))
colnames(genes.df) = "gene"
genes.df$type = NA
genes.df[(grepl("%$", genes.df$gene)),]$type = "wildcard"
genes.df[!(grepl("%$", genes.df$gene)),]$type = "gene"

#create wildcard gene lists
if (dim(genes.df[(genes.df$type == "wildcard"),])[1] > 0 ){
  wildcards = paste("ens.gene_name LIKE '", paste(as.character(genes.df[(genes.df$type == "wildcard"),]$gene), 
                                                  collapse= "' OR ens.gene_name LIKE '"), "'", sep="")
} else
  wildcards = NULL

#create non-wildcard gene list
if (dim(genes.df[(genes.df$type == "gene"),])[1] > 0){
  genes = paste("WHERE (ens.gene_name IN ('", paste(as.character(genes.df[(genes.df$type == "gene"),]$gene), collapse="','")
                , "')", sep="")
} else
  genes = NULL


#if wildcard and genes, add statement to gene list
if ((length(wildcards) > 0) & length(genes) > 0 ){
  gene_list = paste(genes, "OR", wildcards, sep=" ")
} else if 
(length(wildcards) > 0 & length(genes) < 1){
  gene_list = paste("WHERE ", wildcards)
} else if (length(wildcards) < 1 & length(genes) > 0){
  gene_list = genes
} else 
  print("No genes of interest found.")


##################
##  Build query ##
##################
#build illumina query
if (plat_arg == 'illumina' |plat_arg == 'Illumina'){
  ben_query = paste0("SELECT p.sample_id, p.qual, p.filter, k.id as rsID, (k.alle_freq * 100) as kav_pct, k.alle_cnt as kav_count, gene_name, p.chr, p.pos, p.ref, p.alt, p.gt, 
                     (CASE  WHEN SUBSTRING(p.sample_id, -2) = '01' THEN 'M' 
                      WHEN SUBSTRING(p.sample_id, -2) = '02' THEN 'F'   
                      WHEN SUBSTRING(p.sample_id, -2) = '03' THEN 'NB' END) as member, 
                     CONCAT(gene_name, ':', p.chr, ':', CAST(p.pos AS STRING)) as variant_id 
                     FROM 
                     (SELECT DISTINCT p.sample_id, p.qual, p.filter, ens.gene_name, p.chr, p.pos, p.ref, p.alt, p.gt 
                      FROM public_hg19.ensembl_genes ens, ", db_arg, " AS p ", 
                     gene_list, 
                     ") AND ens.chromosome NOT LIKE 'H%' 
                     AND ens.chromosome = p.chr 
                     AND (p.pos >= ens.start AND p.pos <= ens.stop) 
                     AND p.gt IS NOT NULL  ) AS p 
                     LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k 
                     ON p.chr = k.chromosome 
                     AND p.pos = k.pos 
                     AND p.ref = k.ref 
                     AND p.alt = k.alt 
                     WHERE (k.alle_freq < .", kav_arg, " OR k.alle_freq IS NULL)")
  
  query = strwrap(ben_query, width=nchar(ben_query), simplify=TRUE)
}

#build comgen query
if (plat_arg == 'comgen' |plat_arg == 'cgi' | plat_arg == 'Comgen' |plat_arg == 'Cgi'){
  ben_query = paste0("SELECT p.sample_id, p.allele1varquality, p.totalreadcount, k.id as rsID, (k.alle_freq * 100) as kav_pct, 
k.alle_cnt as kav_count, gene_name, p.chr, p.start, p.stop, p.ref, p.allele1seq, p.allele2seq, p.zygosity,  
    (CASE 
    WHEN SUBSTRING(p.sample_id, -1) = 'M' THEN 'M' 
    WHEN SUBSTRING(p.sample_id, -1) = 'F' THEN 'F' 
    WHEN SUBSTRING(p.sample_id, -2) = 'NB' THEN 'NB'
    END) as member, CONCAT(gene_name, ':', p.chr, ':', CAST(p.start AS STRING), ':', CAST(p.stop AS STRING)) as variant_id  
FROM 
    (SELECT DISTINCT p.sample_id, p.allele1varquality, p.totalreadcount, ens.gene_name, p.chr, p.start, p.stop, p.ref, p.allele1seq, p.allele2seq, p.zygosity
    FROM public_hg19.ensembl_genes ens, ", db_arg, " AS p ", 
                     gene_list, 
                     ") AND ens.chromosome NOT LIKE 'H%' 
    AND ens.chromosome = p.chr 
    AND (p.start >= ens.start AND p.stop <= ens.stop)
    AND p.zygosity IS NOT NULL) AS p
LEFT JOIN /* +SHUFFLE */ public_hg19.kaviar k
                     ON p.chr = k.chromosome
                     AND (k.pos BETWEEN p.start and p.stop)
                     AND p.ref = k.ref
                     AND p.allele1seq = k.alt
                     WHERE (k.alle_freq < .", kav_arg, " OR k.alle_freq IS NULL)")

query = strwrap(ben_query, width=nchar(ben_query), simplify=TRUE)
}

#########################
## Run Query on Impala ##
#########################
query_results = sqlQuery(conn, query)

#################################
## structure data for analysis ##
#################################
#split apart data frame by member id and examine data
mom = query_results[which(query_results$member == "M"),]
dad = query_results[which(query_results$member == "F"),]
nb = query_results[which(query_results$member == "NB"),]

#######################
#find homozygous alt ##
#######################
#merge together to find genes in parents that match nb variants on each gene, by chr and position
nb_mom_hom_alt = merge(hom_alt_nb, het_mom, by="variant_id")

#merge nb with dad
nb_dad_hom_alt= merge(nb_mom_hom_alt, het_dad, by="variant_id")

#merge mom and dad's matching varaints
hom_alts = merge(nb_mom_hom_alt, nb_dad_hom_alt, by="variant_id")

#get rid of unnecessary columns
hom_alt.df = hom_alts[,c(1:13,25:26,38:39)]
colnames(hom_alt.df) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_freqPct", 
                         "kav_count", "gene_name", "chr", "pos","ref", "nb_alt", "nb_gt", 
                         "m_alt", "m_gt", "f_alt","f_gt")

#find variants where mother, father and newborn have matching alts at same chr/pos, or
#mom and dad alts match and nb call is null
hom_alt.df = hom_alt.df[which((hom_alt.df$nb_alt == hom_alt.df$m_alt) & (hom_alt.df$nb_alt == hom_alt.df$f_alt) | 
                                (hom_alt.df$m_alt == hom_alt.df$f_alt & hom_alt.df$nb_alt == "NULL")
),]

#############################
#find homozygous reference ##
#############################
hom_ref_nb = nb[which(nb$gt == "0/0"),]

#merge hom ref nb with 19 het mom variants, matching by chr and position
hom_ref_nb_mom = merge(hom_ref_nb, het_mom, by="variant_id")

#merge hom ref nb with het dad
hom_ref_nb_dad = merge(hom_ref_nb, het_dad, by="variant_id")

#merge mom and dad matching variants
hom_refs = merge(hom_ref_nb_mom, hom_ref_nb_dad, by= "variant_id")

#get rid of unnecessary columns
hom_ref.df = hom_refs[,c(1:13,25:26,38:39)]
colnames(hom_ref.df) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_freqPct",
                         "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt", "nb_gt",
                         "m_alt", "m_gt", "f_alt","f_gt")

#find variants where mother, father and newborn have matching alts at same chr/pos, or
#mom and dad alts match and nb call is null
hom_ref.df = hom_ref.df[which((hom_ref.df$nb_alt == hom_ref.df$m_alt) & (hom_ref.df$nb_alt == hom_ref.df$f_alt) | 
                                (hom_ref.df$m_alt == hom_ref.df$f_alt & hom_ref.df$nb_alt == "NULL")
),]

#####################
#find compound het ##
#####################
#subset for possible zygosity for compound het
het_nb = nb[which(nb$gt == "0/1"),]

#function to check for comp_hets
find_comphet = function(x){
  #check if more than one variant found per gene
  if (length(unique(x$pos)) >1 ) {
    #find matching variants from each parent
    mom_vars = mom[grep(paste(x$variant_id, collapse="|"), mom$variant_id),]
    dad_vars = dad[grep(paste(x$variant_id, collapse="|"), dad$variant_id),]
    #if either parent has a variant the other does not 
    #and there is at least one variant from each, mark as comp_het
    if ((dim(mom_vars[!(mom_vars$variant_id %in% dad_vars$variant_id),])[1] > 0 | 
           dim(dad_vars[!(dad_vars$variant_id %in% mom_vars$variant_id),])[1] > 0) &
          (dim(mom_vars)[1] > 0 & dim(dad_vars)[1] > 0)){
      rbind(mom_vars, dad_vars)
    } else 
      print("Variants at different positions on the gene are from the same parent.")     
  }else #(length(unique(x$pos)) >1 ) = false
    print("Variants on gene only found in same position.")
}#end of function

##apply function to newborn comp het candidates
comp_hets = by(het_nb, as.character(het_nb$gene_name), find_comphet)

#data frame results - automate this dammit Summer
comp_het.df = rbind(comp_hets$HSD17B1, comp_hets$HSD17B13, comp_hets$HSD17B4, comp_hets$HSD17B7P2)
#combine with newborn variants
comp_het.df = merge(het_nb, comp_het.df, by="variant_id")
comp_het.df = comp_het.df[,c(1:13,25:27)]
colnames(comp_het.df) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_pct",
                          "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt",
                          "nb_gt", "parent_alt", "parent_gt", "parent")

########################
## Find MIE in Het NB ##
#########################
nb_het = nb[which(nb$gt == "0/1"),]

#find matching variants from each parent
mom_vars = mom[grep(paste(nb_het$variant_id, collapse="|"), mom$variant_id),] #10 vars
dad_vars = dad[grep(paste(nb_het$variant_id, collapse="|"), dad$variant_id),] #9 vars

#merge together to find intersection of gene:chr:pos
mie_het_nb_mom = merge(nb_het, mom_vars, by = "variant_id")
mie_hets_nb_dad = merge(nb_het, dad_vars, by = "variant_id")
mie_hets = merge(mie_het_nb_mom, mie_hets_nb_dad, by= "variant_id")

#get rid of unnecessary columns
mie_hets = mie_hets[,c(1:13,25:26,38:39)]
colnames(mie_hets) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_freqPct",
                       "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt", "nb_gt",
                       "m_alt", "m_gt", "f_alt","f_gt")

#variants not in mie_hets and not marked as comp_het are missing info from one parent
missing_info = nb_het[grep(paste(mie_hets$variant_id, collapse="|"), nb_het$variant_id, invert=TRUE),]
missing_info = missing_info[grep(paste(comp_het.df$variant_id, collapse="|"), missing_info$variant_id, invert=TRUE),]

#find variants that are in congruence with MI laws
no_mie = rbind(
  #if the mother is 0/0 and father is 0/1, then the nb alt must match the father's alt or be null
  mie_hets[(mie_hets$m_gt == "0/0" & mie_hets$f_gt == "0/1" & (mie_hets$nb_alt == mie_hets$f_alt | mie_hets$nb_alt == "NULL")),],
  #if the father is 0/0 and the mother is 0/1, then the nb alt must match the mother's alt or be null
  mie_hets[(mie_hets$f_gt == "0/0" & mie_hets$m_gt == "0/1" & (mie_hets$nb_alt == mie_hets$m_alt| mie_hets$nb_alt == "NULL")),],
  #if both parents are het, nb alt and parent alts must be the same, or nb alt is null
  mie_hets[((mie_hets$m_gt == "0/1" & mie_hets$f_gt== "0/1") & (mie_hets$nb_alt == mie_hets$f_alt & mie_hets$nb_alt == mie_hets$m_alt)| (mie_hets$m_alt == mie_hets$f_alt & mie_hets$nb_alt == "NULL")),],
  #if the mother is 1/1, then the father must be 0/1 or 0/0 and the nb alt must match the mother's alt
  #or be null
  mie_hets[(mie_hets$m_gt == "1/1" & (mie_hets$f_gt == "0/1" | mie_hets$f_gt == "0/0") & (mie_hets$nb_alt == mie_hets$m_alt & mie_hets$nb_alt == mie_hets$f_alt | mie_hets$m_alt == mie_hets$f_alt & mie_hets$nb_alt == "NULL")),],
  #if the father is 1/1, then the mother must be 0/1 or 0/0 and the nb alt must match the fathers's alt
  mie_hets[(mie_hets$f_gt == "1/1" & (mie_hets$m_gt == "0/1" | mie_hets$m_gt == "0/0") & (mie_hets$nb_alt == mie_hets$f_alt & mie_hets$nb_alt == mie_hets$m_alt | mie_hets$f_alt == mie_hets$m_alt & mie_hets$nb_alt == "NULL")),])

#if mie_het variants are not in no_mie set, then they are potentially MIE
het_mie = mie_hets[grep(paste(no_mie$variant_id, collapse="|"), mie_hets$variant_id, invert=TRUE),]

###################################
## Find MIE in homozygous ref nb ##
###################################
nb_hom_ref = nb[which(nb$gt == "0/0"),]

#find equivalent variants in parents
mom_hom_vars = na.omit(mom[match(nb_hom_ref$variant_id, mom$variant_id),])
dad_hom_vars = na.omit(dad[match(nb_hom_ref$variant_id, dad$variant_id),])

#merge together to find matching parent varients
mie_hom_nb_mom = merge(nb_hom_ref, mom_hom_vars, by = "variant_id")
mie_hom_nb_dad = merge(mie_hom_nb_mom, dad_hom_vars, by = "variant_id")
mie_homs = merge(mie_hom_nb_mom, mie_hom_nb_dad, by= "variant_id")

#clean up resuts
mie_homs = mie_homs[,c(1:13,25:26,38:39)]
colnames(mie_homs) = c("variant_id", "sample_id", "qual", "filter", "rsID", "kav_freqPct",
                       "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt", "nb_gt",
                       "m_alt", "m_gt", "f_alt","f_gt")

#variants that are not in mie_homs and not marked hom ref are missing parent info
missing_info_miehoms = nb_hom_ref[grep(paste(mie_homs$variant_id, collapse="|"), nb_hom_ref$variant_id, invert=TRUE),]
missing_info_miehoms = missing_info[grep(paste(hom_ref.df$variant_id, collapse="|"), missing_info$variant_id, invert=TRUE),]
missing_info_miehoms

#find variants that are in congruence with MI laws
no_mie_homs = rbind(
  #if the mother is 0/1, then the dad must be 0/0 
  mie_homs[(mie_homs$m_gt == "0/1" & mie_homs$f_gt == "0/0"),],
  #if the father is 0/1, then the mother must be 0/0
  mie_homs[(mie_homs$d_gt == "0/1" & mie_homs$m_gt == "0/0"),],
  #both parents are 0/0
  mie_homs[(mie_homs$m_gt == "0/0" & mie_homs$f_gt== "0/0"),]
)

#if nb_hom variants are not in no_mie_homs set, then they are MIE
hom_ref_mie = na.omit(nb_hom_ref[!match(nb_hom_ref$variant_id, no_mie_homs$variant_id),])

###################################
## Find MIE in homozygous alt nb ##
###################################
nb_hom_alt = nb[which(nb$gt == "1/1"),]

#find equivalent variants in parents
mom_hom_alt = na.omit(mom[match(nb_hom_alt$variant_id, mom$variant_id),])
dad_hom_alt = na.omit(dad[match(nb_hom_alt$variant_id, dad$variant_id),])

#merge together to find intersection of gene:chr:pos
mie_alt = merge(nb_hom_alt, mom_hom_alt, by = "variant_id")
mie_alt = merge(mie_alt, dad_hom_alt, by = "variant_id")

#clean up resuts
mie_alts = mie_alt[,c(1:13,25:26,38:39)]
colnames(mie_alts) = c("variant_id", "qual", "filter", "rsID", "kav_pct",
                       "kav_count", "gene_name", "chr", "pos", "ref", "nb_alt",
                       "nb_gt", "m_alt", "m_gt", "f_alt", "f_gt")

#find variants that are in congruence with MI laws
mie_hom_alt = rbind(
  #if the mother is 0/1, then the father must be 1/1
  mie_alts[(mie_alts$m_gt == "0/1" & mie_alts$f_gt == "1/1"),],
  #if the father is 0/1, then the mother must be 1/1
  mie_alts[(mie_alts$d_gt == "0/1" & mie_alts$m_gt == "1/1"),],
  #both parents are 1/1
  mie_alts[(mie_alts$m_gt == "1/1" & mie_alts$f_gt == "1/1"),]
)

#if nb_hom variants are not in no_mie_homs set, then they are MIE
hom_alt_mie = na.omit(nb_hom_alt[!match(nb_hom_alt$variant_id, mie_hom_alt$variant_id),])

####################################
## Merge Results and save to file ##
####################################

#remove variant id

#comp hets, find both parents alts for each position


comp_het_nb.df = comp_het.df[,c(1:13)]


mom_comp_hets = na.omit(mom[match(comp_het_nb.df$variant_id, mom$variant_id),])
dad_comp_hets = na.omit(dad[match(comp_het_nb.df$variant_id, dad$variant_id),])

nb_mom_comp_hets = merge(comp_het_nb.df, mom_comp_hets, by="variant_id")
nb_mom_comp_hets = nb_mom_comp_hets[,c(1:13, 24:25)]

nb_dad_comp_hets = merge(comp_het_nb.df, nb_mom_comp_hets, by="variant_id")
nb_dad_comp_hets = nb_dad_comp_hets[,c(1:13, 24:25)]
head(nb_dad_comp_hets)

comp_hets.df = merge(nb_mom_comp_hets, nb_dad_comp_hets, by="variant_id")

colnames(comp_hets.df)

test = comp_hets[,c(2:13,24:25,)]
head(test)

head(hom_alt.df)
head(comp_het.df)




