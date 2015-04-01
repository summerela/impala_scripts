library(VariantAnnotation)

#specify location of vcf file
vcf_in = "~/Downloads/clinvar.vcf.gz"
#specify location to write info field file
info_out = "~/impala_scripts/"
                                           
#################################################
## Grab info fields and description for db def ##
#################################################
get_info = function(vcf){
  #examine vcf header to grab info fields
  hdr =  scanVcfHeader(vcf)
  #build data frame from info fields
  info = data.frame(field = rownames(info(hdr)), type="STRING", comment="COMMENT", description = paste("'", info(hdr)$Description, "'", sep=""))
  ##write info to a file for db_def
  write.table(info, paste(info_out, paste(lapply(strsplit(basename(vcf), "[.]"), function(x) x[1]), "info", sep="_"), ".txt", sep=""), 
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}

get_info(vcf_in)

######################################
## check for fields with delimiters ##
######################################
#read in only info fields
param = ScanVcfParam(fixed=NA, samples=NA)

#read in vcf file with paramaters
vcf_data = readVcf(vcf_in, "hg19", param=param)

#turn info field into data frame for searching 
vcf_info = as.data.frame(info(vcf_data))

apply(as.matrix(vcf_info), 1, function(x){ sum(grepl("|", x, perl=T)) > 0 })

test_col = "CLNSIG"

if (class(test_col) = 'character')
  {
  paste("yep its true")
}



