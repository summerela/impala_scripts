#!/usr/bin/env Rscript

#for windows, export R path
system(paste("export PATH=$PATH:/usr/bin/R"), wait=TRUE)

###############################
## Read in command line args ##
###############################

##Grab user args from command line
args=(commandArgs(TRUE))

##print help if both args not provided 
if(length(args) < 2) {
  cat("\n **Enter arguments for both gene and db:** \n")
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      variant_analyze.R annotate variants then test for zygosity and MIE
      
      Arguments:
      --gene  - Genes of interest, surrounded by quotes, comma-sep, no spaces
      --db      - Database to search in, ex.p7_ptb.illumina_variant
      --help             - print help text
      
      Example:
      ./variant_analyze.R --gene='BRCA1' --db='p7_ptb.illumina_variant' \n\n")
 
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




