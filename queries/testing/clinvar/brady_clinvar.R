brady = read.csv("./brady_clinvar.txt", sep="\t", header=TRUE)

head(brady)

brady$chrom = gsub("chr", "", lapply(strsplit(as.character(brady$position), ":"), function(x) x[1]))
brady$trio_id = lapply(strsplit(as.character(brady$identifier.or.consent), "-"), function(x) x[3])

trios = c("M", "F", "NB")
brady[!which(brady$trio_id %in% trios),]

#homozygous
#chr8
#clinsig = 5
#mother or father

sigs = brady[which(brady$chrom == "8" & brady$zygosity == "homozygous" 
                   & brady$trio_id != "NB" & brady$clinvar.pathogenicity == "5"),]

head(sigs)

