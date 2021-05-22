## author : ambati@stanford.edu
## convert DR1 alleles to amino acid calls
arg=commandArgs(trailingOnly = T)
#arg[1] = 'DRB1'
#arg[2] = 'outFile'
PackList = rownames(installed.packages())
if('data.table' %in%  PackList){
  require(data.table)
} else {
  install.packages('data.table')
  require(data.table)
  
}

if('HIBAG' %in% PackList){
  require(HIBAG)
} else {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("HIBAG")
  require(HIBAG)
}

## main function 
HLA2AA=function(arg){
  hlaFile=fread(arg[1])
  hlaObj=hlaAllele(sample.id = hlaFile$sample.id, H1 = hlaFile$DRB1.1, H2=hlaFile$DRB1.2, prob = hlaFile$DRB1.prob, locus = 'DRB1')
  hla.aa=hlaConvSequence(hla = hlaObj, code = "P.code.merge")
  filtered_HLA_gene <- hla.aa$value[hla.aa$value$prob > 0.3,]
  pos.table = summary(hla.aa)
  increment <- pos.table[1,"Pos"] - 1
  pos.table[,"Pos"] <- pos.table[,"Pos"] - increment
  hla_pos <-  pos.table[,"Pos"]
  out=lapply(hla_pos, function(pos){
    HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
    HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("ID","allele")
    names(HLA_allele2) <- c("ID","allele")
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)
    HLA_aa_code <- HLA_allele
    HLA_aa_code$allele <- substr(HLA_aa_code$allele, pos, pos)
    HLA_count <- as.data.frame.matrix(table(HLA_aa_code))
    names(HLA_count)[names(HLA_count) == "-"] <- "Ref"
    names(HLA_count)[names(HLA_count) == "*"] <- "Amb"	
    names(HLA_count)[names(HLA_count) == "."] <- "CNV" 
    names(HLA_count)=paste0(names(HLA_count), '_', pos+increment)
    head(HLA_count)
  })
  return(out)
}
## convert to HLA amino acids from HLA -DR calls
convList=HLA2AA(arg)
convDF = data.table(sample.id=rownames(convList[[1]]))
for(pos in convList){
  temp=data.table(pos, keep.rownames = T)
  if(identical(temp$rn, convDF$sample.id)){
    convDF = cbind.data.frame(convDF, pos)
  }
}
setDT(convDF)
convDF
## write out the amino acid DF
fwrite(convDF, file = arg[2])
