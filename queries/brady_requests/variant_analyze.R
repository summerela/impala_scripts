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

args = "HOX%,BRCA1,BRCA2,RAS%"

args.df = as.data.frame(matrix(unlist(strsplit(args, ","))))
colnames(args.df) = "gene"
args.df$type = NA
args.df[(grepl("%$", args.df$gene)),]$type = "wildcard"
args.df[!(grepl("%$", args.df$gene)),]$type = "gene"


paste("ens.gene_name LIKE '", paste(as.character(args.df[(args.df$type == "wildcard"),]$gene), 
                                  collapse= "' OR ens.gene_name LIKE '"), "'", sep="")

paste("WHERE ens.gene_name IN ('", paste(as.character(args.df[(args.df$type == "gene"),]$gene), collapse="','")
                                        , "')", sep="")




#######################
##  Run impala query ##
#######################
query = paste("SELECT * FROM public_hg19.ensembl_genes ens WHERE ens.gene_name = '", argsL$gene, "' LIMIT 5", sep="")

print(sqlQuery(conn, query))

