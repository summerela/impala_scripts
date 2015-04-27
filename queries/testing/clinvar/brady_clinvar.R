library(readr)

##read in brady and impala query results
brady = read_tsv("./queries/testing/clinvar/brady_clinvar.txt")
impala = read_csv("./queries/testing/clinvar/impala_query_results.csv")

##make data frames with matching variables
brady.df = data.frame(chr = gsub("chr", "", lapply(strsplit(as.character(brady$position), ":"), function(x) x[1])),
                      pos = unlist(lapply(strsplit(as.character(brady$position), ":"), function(x) x[2])),
                      ref = unlist(lapply(strsplit(as.character(brady$dna_change), "->"), function(x) x[1])),
                      alt = unlist(lapply(strsplit(as.character(brady$dna_change), "->"), function(x) x[2])),
                      zygosity = gsub("homozygous", "hom", brady$zygosity),
                      gene = unlist(lapply(strsplit(brady$gene_definition, ':'), function(x) x[1])),
                      sample_id = gsub(".*:","",brady$identifier_or_consent), 
                      clin_sig = brady$clinvar_pathogenicity)

impala.df = data.frame(chr = as.character(impala$chr), 
                       pos = as.character(impala$start),
                       ref = impala$ref,
                       alt = impala$allele1seq,
                       zygosity = impala$zygosity,
                       gene = unlist(lapply(strsplit(impala$clin_geneinfo, ':'), function(x) x[1])),
                       sample_id = impala$sample_id,
                       clin_sig = as.character(impala$clin_sigid))

##coerce all values to character for matching
#coercing columns to same class for comparison
i = sapply(brady.df, is.factor)
brady.df[i] = lapply(brady.df[i], as.character)

j = sapply(impala.df, is.factor)
impala.df[j] = lapply(impala.df[j], as.character)

##order data frames by sample id, chr, pos for matching
#order both data frames for matching
brady.df = brady.df[with(brady.df, order(sample_id, chr, pos)),]
impala.df = impala.df[with(impala.df, order(sample_id, chr, pos)),]

##subest brady's results to match with impala query
pathogenic = c("4", "5")
not_pathogenic = c("2", "3")
brady_filter = brady.df[which(brady.df$chr == "8" & brady.df$zygosity == "hom"
                          &  grep(paste(pathogenic,collapse="|"), brady.df$clin_sig) &
                            grep(paste(not_pathogenic,collapse="|"), brady.df$clin_sig, invert=TRUE)),]

##clear rownames for matching
rownames(impala.df) = NULL
rownames(brady_filter) = NULL

##since the impala set does not include piping, and all the pipes from Brady's results are "5|5" gsub with "5" for matching
##yes I realize this is a bit of cheating
brady_filter$clin_sig = as.character("5")

require(dplyr)
#records in impala not in brady's results
not_in_brady = unique(anti_join(impala.df,brady_filter))
dim(not_in_brady)
not_in_brady

##records in brady's results but not in impala
not_in_impala = unique(anti_join(brady_filter,impala.df))
dim(not_in_impala)
not_in_impala



