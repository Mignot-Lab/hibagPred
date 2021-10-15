## read HLA
fileIO=function(filePre){
  pathTohla=list.files(path = 'hlaOut/', pattern = paste0(filePre, '*'), full.names = T)
  if (length(pathTohla) >= 3){
    message(paste0('IDENTIFIED HLA IMPUTE FILES IN hlaOut/ ', pathTohla, ' ', timestamp(), collapse = '\n'), appendLF = T)
    hlaNames = gsub('IMP_|.csv', '', unique(str_extract(pattern = "IMP_[^\\s]+", pathTohla)))
    hlaDFlist = lapply(pathTohla, function(hlaFile){
      temp = fread(hlaFile, key='sample.id', select = c(1:4)) # take first 4 cols to output
    })
    names(hlaDFlist) = hlaNames
    return(hlaDFlist)
  } else {
    message('HAVE IMPUTATIONS BEEN PERFORMED ? RUN HibagPred.R FIRST OR CHECK DIR LOCATION hlaOut/ ', timestamp() )
    stop()
  }
}

## amino acid conv 
HLA2AA=function(arg){
  hlaFile=fread(arg[1])
  locus = arg[3]
  hlaObj=hlaAllele(sample.id = hlaFile[[1]], H1 = hlaFile[[2]], H2=hlaFile[[3]], prob = hlaFile[[4]], locus = locus, max.resolution = '4-digit')
  hla.aa=hlaConvSequence(hla = hlaObj, code = "P.code.merge")
  filtered_HLA_gene = hla.aa$value[hla.aa$value$prob > 0.3,] # place holder
  pos.table = summary(hla.aa)
  increment = pos.table[1,"Pos"] - 1
  pos.table[,"Pos"] = pos.table[,"Pos"] - increment
  hla_pos =  pos.table[,"Pos"]
  out=lapply(hla_pos, function(pos){
    HLA_allele1 = filtered_HLA_gene[, c("sample.id","allele1")]
    HLA_allele2 = filtered_HLA_gene[, c("sample.id","allele2")]
    names(HLA_allele1) = c("ID","allele")
    names(HLA_allele2) = c("ID","allele")
    HLA_allele = rbind(HLA_allele1, HLA_allele2)
    HLA_aa_code = HLA_allele
    HLA_aa_code$allele = substr(HLA_aa_code$allele, pos, pos)
    HLA_count = as.data.frame.matrix(table(HLA_aa_code))
    names(HLA_count)[names(HLA_count) == "-"] = "Ref"
    names(HLA_count)[names(HLA_count) == "*"] = "Amb"	
    names(HLA_count)[names(HLA_count) == "."] = "CNV" 
    names(HLA_count)=paste0(locus, '_', names(HLA_count), '_', pos+increment)
    #HLA_count$PROB = hlaFile[[4]]
    return(HLA_count)
  })
}

##function to make haplotypes class II
makeHaps=function(hlaDFlist){
  outPutList = NULL
  hlasNeeded = c('DQB1', 'DQA1', 'DRB1')
  if(all(hlasNeeded %in% names(hlaDFlist))){
    dr1.df = hlaDFlist[['DRB1']][, 1:3]
    dqa.df = hlaDFlist[['DQB1']][, 1:3]
    dqb.df = hlaDFlist[['DQA1']][, 1:3]
    setkey(dr1.df, key='sample.id');setkey(dqb.df, key='sample.id');setkey(dqa.df, key='sample.id')
    genos=as.data.frame(dqa.df[dqb.df][dr1.df])#
    if(dim(genos)[1] > 1){
    #locus.label = c('DQA1', 'DQB1', 'DRB1')
      outPutList$IDS = genos$sample.id
      fit=haplo.em(geno = genos[, -1], locus.label=c('DQB1', 'DQA1', 'DRB1'))
      outPutList$hapFreq = print(fit)
      outPutList$genos = genos
      haplCalls = data.table(sample.id=genos$sample.id[fit$subj.id], hap.1=fit$hap1code,  hap.2=fit$hap2code, hap.prob=fit$post)
      haplotypes = apply(fit$haplotype, 1, function(x) paste0(x, collapse = "_"))
      ##assign actuall lhaps
      haplCalls$hap.1=haplotypes[haplCalls$hap.1]
      haplCalls$hap.2=haplotypes[haplCalls$hap.2]
      outPutList$fit = fit
      outPutList$haplCalls = haplCalls#hlaAllele(sample.id = haplCalls$id, H1 = haplCalls$hap1, H2=haplCalls$hap2, prob = haplCalls$prob)
      #outPutList$haplCalls$locus = 'DQB_DQA_DRB'
      outPutList$haplotypes = haplotypes
      return(outPutList)
    } else {
      message('CHECK OBJECT FOR PRESENCE OF DRB1 DQA1 & DQB1 HLA CALLS', timestamp())
      return(NULL)
    }
  }
}

## function to parse PCS and DX 
metaIO=function(metaFile){
  metaIn = fread(metaFile, key='sample.id')
  if ('Pheno' %in% names(metaIn)){
    message(paste0('PHENO NAME IDENTIFIED AS PHENO', timestamp()))
    phenoVec = as.numeric(metaIn[['Pheno']])
    casesN = sum(phenoVec == 1)
    controlsN = sum(phenoVec == 0)
    if (casesN > 0 & controlsN > 0){
      message(paste0('IDENTIFIED CASES N = ', casesN))
      message(paste0('IDENTIFIED CTRLS N = ', controlsN))
    } else {
      message('HAVE CASES BEEN CODED AS 1 AND CONTROLS AS 0 ? CODING PLINK PHENO TO 0- CONTROLS AND 1 - CASES', timestamp())
      phenoVec = ifelse(phenoVec == 1, 0, ifelse(phenoVec == 2, 1, NA))
      casesN = sum(phenoVec == 1)
      controlsN = sum(phenoVec == 0)
      message(paste0('IDENTIFIED CASES N = ', casesN))
      message(paste0('IDENTIFIED CTRLS N = ', controlsN))
      metaIn[['Pheno']] = phenoVec
    }   
    return(metaIn)
  } else {
    stop('META FILE DOESNT CONTAIN PHENO NAME INPUT')
  }

}
