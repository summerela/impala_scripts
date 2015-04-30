#load libraries with no whining
suppressMessages(library(readr))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(AnnotationHub))
suppressMessages(library(biomaRt))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(RODBC))

#####################################################
## connect to impala and pull down illumina data  ### 
#####################################################

#remove extraneous columns
ill = read_csv("~/GitHub/impala_scripts/queries/brady_requests/itmi_puzzle_data.csv", 
               col_types = "ci_cci________________________________c__c")

#add trio and member id columns
ill$trio_id = unlist(lapply(strsplit(as.character(ill$sample_id), "-"), function(x) paste(x[1], x[2], sep="-")))
ill$member_id =  unlist(lapply(strsplit(as.character(ill$sample_id), "-"), function(x) x[3]))

#subset for het mothers, het fathers and all newborns
subset = ill[which(ill$member_id == "03" | ill$member_id == "02" & ill$gt == "0/1"| data$member_id == "01" & data$gt == "0/1"),]

##############################
## find ensembl annotations ##
##############################

annotate_vars  = function(variant_df, gene_list){
  #load hg19 annotation data
  hub <- AnnotationHub()
  hum = hub[['AH10684']]
  #load ensembl data from biomart
  ensembl = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org")
  #create granges object from subset
  variant_ranges = with(variant_df, GRanges(seqnames=chrom, 
                                        ranges=IRanges(pos, width=1, 
                                        names=paste(sample_id, ":",ref, ">", alt, sep="")),
                                        sample_id = sample_id,
                                        ref = ref,
                                        alt = alt,
                                        qual = qual,
                                        gt = dpf,
                                        trio_id = trio_id,
                                        member_id = member_id))
  #map variants to ensembl
  loc <- locateVariants(query=variant_ranges, subject=ensembl, region=AllVariants())
  #exmaine regions returned
  table(loc$LOCATION)
  #data frame results
  loc.df = as.data.frame(loc, row.names= seq(1:length(loc)))
  #merge results with variant subset
  annots.df = merge(loc.df, variant_df, by.x=c("seqnames", "start"), by.y=c("chrom", "pos"))
  #filter out intergenic regions
  annots.df = annots.df[(annots.df$LOCATION != "intergenic"),]
  #map ensembl geneid to gene name
  cols <- c("SYMBOL", "UNIPROT", "REFSEQ", "GENENAME")
  keys <- unique(annots.df$GENEID)
  gene_names = select(org.Hs.eg.db, keys, cols, keytype="ENSEMBL")
  #final list of annotated variants
  variant_annot = merge(annots.df, gene_names, by.x="GENEID", by.y="ENSEMBL")
  candidates = variant_annot[which(gene_list %in% variant_annot$SYMBOL),]
  candidates
}

#list genes of interest
genes_of_interest =  c('HADH', 'HADHA', 'HADHB', 'ACAA1', 'ACAA2', 'EHHADH', 'ECHS1')

#annotate
annotate_vars(subset, genes_of_interest)







