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
HLA2AA=function(hlaFile, locus, probCutoff){
  hlaObj=hlaAllele(sample.id = hlaFile[[1]], H1 = hlaFile[[2]], H2=hlaFile[[3]], prob = hlaFile[[4]], locus = locus, max.resolution = '4-digit')
  hla.aa=hlaConvSequence(hla = hlaObj, code = "P.code.merge")
  filtered_HLA_gene = hla.aa$value[hla.aa$value$prob > probCutoff,] # place holder
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
      haplCalls=haplCalls[, .SD[which.max(hap.prob)], by=sample.id]
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

aminoDRB1 = function(hlaFile, locus, probCutoff){
  hlaObj=hlaAllele(sample.id = hlaFile[[1]], H1 = hlaFile[[2]], H2=hlaFile[[3]], prob = hlaFile[[4]], locus = locus, max.resolution = '4-digit')
  hla.aa=hlaConvSequence(hla = hlaObj, code = "P.code.merge")
  filtered_HLA_gene = hla.aa$value[hla.aa$value$prob > probCutoff,]
  pos.table = summary(hla.aa)
  increment = pos.table[1,"Pos"] - 1
  pos.table[,"Pos"] = pos.table[,"Pos"] - increment
  hla_pos =  pos.table[,"Pos"]
  posNeeded = c(13, 33, 57)-increment
  lapply(posNeeded, function(pos){
    a1 = substr(filtered_HLA_gene$allele1, pos, pos)
    a2 = substr(filtered_HLA_gene$allele2, pos, pos)
    a1[a1 == "-"] = "Ref";a1[a1 == "*"] = "Amb";a1[a1 == "."] = "CNV"
    a2[a2 == "-"] = "Ref";a2[a2 == "*"] = "Amb";a2[a2 == "."] = "CNV"
    alleleNames = paste0(locus, '_', pos+increment)
    a1 = paste0(pos+increment, '_', a1);a2 = paste0(pos+increment, '_', a2)
    outTemp=data.table(sample.id = filtered_HLA_gene$sample.id, a1, a2)
    names(outTemp)[2:3] = paste0(alleleNames, '.', 1:2)
    return(outTemp)
  })
}

## Function to construct AminoAcid Haplotypes
aminoHap = function(DRB1AAList){
    genos = data.frame(sample.id = DRB1AAList[[1]]$sample.id)
    for(df in DRB1AAList){
      if(identical(df$sample.id, genos$sample.id)){
        genos=cbind.data.frame(genos, df[, 2:3])
      } else {
        stop('SAMPLE IDS DONT MATCH')
      }
    }
    if(dim(genos)[1] > 20){
      timestamp()
      message('FITTING HAPLOTYPES TO DRB1 AMINOACID HAPLOTYPES')
      outPutList = NULL
      outPutList$IDS = genos$sample.id
      fit=haplo.em(geno = genos[, -1], locus.label=c('DRB1_13', 'DRB1_33', 'DRB1_57'))
      outPutList$hapFreq = print(fit)
      outPutList$genos = genos
      haplCalls = data.table(sample.id=genos$sample.id[fit$subj.id], hap.1=fit$hap1code,  hap.2=fit$hap2code, hap.prob=fit$post)
      haplCalls=haplCalls[, .SD[which.max(hap.prob)], by=sample.id]
      haplotypes = apply(fit$haplotype, 1, function(x) paste0(x, collapse = "_"))
      ##assign actuall lhaps
      haplCalls$hap.1=haplotypes[haplCalls$hap.1]
      haplCalls$hap.2=haplotypes[haplCalls$hap.2]
      outPutList$fit = fit
      outPutList$haplCalls = haplCalls#hlaAllele(sample.id = haplCalls$id, H1 = haplCalls$hap1, H2=haplCalls$hap2, prob = haplCalls$prob)
      #outPutList$haplCalls$locus = 'DQB_DQA_DRB'
      outPutList$haplotypes = haplotypes
      #one hot encode the haplotypes
      outPutList$onc = as.data.frame.matrix(table(data.table(sample.id=rep(haplCalls$sample.id, 2), allele = c(haplCalls$hap.1, haplCalls$hap.2)))) %>% data.table(keep.rownames = T)
      return(outPutList$onc)
    } else {
      message('MIN SAMPLE SIZE IS 20')
      return(NULL)
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

allele.Glm.Fishers = function(hla, hla.meta, HLA_count){
  HLA_allele_name = names(HLA_count)[-1]
  if (length(HLA_allele_name) >= 1){ ## more than one allele atleast
    setDT(HLA_count, key='sample.id', keep.rownames = T)
    HLA_count[hla.meta, c('Pheno', 'PC1', 'PC2', 'PC3', 'PC4'):=list(Pheno, PC1, PC2, PC3, PC4), by=.EACHI]
    HLA_count = HLA_count[!is.na(Pheno)]
    n.cases = nrow(HLA_count[Pheno==1])#$sample.id
    n.controls = nrow(HLA_count[Pheno==0])
    HLA_countCarrier = copy(HLA_count)
    HLA_countCarrier[, HLA_allele_name] = data.table(apply(HLA_countCarrier[, HLA_allele_name, with =F], 2, function(x) ifelse(x > 0, 1, 0)))
    lapply(HLA_allele_name, function(allele){
      model.design=paste0("Pheno ~ ",allele,"+PC1+PC2+PC3+PC4")
      fit = glm(model.design, data = HLA_countCarrier, family = 'binomial')
      fitDF = data.table(tidy(fit))
      fitDF=fitDF[term == allele]
      #print(fitDF)
      Case.Car= apply(HLA_count[Pheno == 1, allele, with =F], 2, FUN = function(x) sum(ifelse(x > 0, 1, 0)))
      Case.Not = n.cases-Case.Car
      Control.Car = apply(HLA_count[Pheno == 0, allele, with =F], 2, FUN = function(x) sum(ifelse(x > 0, 1, 0)))
      Control.Not = n.controls-Control.Car
      totalCarrier = apply(HLA_count[, allele, with =F], 2, FUN = function(x) sum(ifelse(x > 0, 1, 0)))
      ## allele freq
      Case.Alleles= sum(HLA_count[Pheno == 1, allele, with =F])
      Case.Alelles.Not = (2*n.cases)-Case.Alleles
      Control.Alleles= sum(HLA_count[Pheno == 0, allele, with =F])
      Control.Alelles.Not = (2*n.controls)-Control.Alleles
      Mat.Carrier = matrix(c(Case.Car,  Control.Car, Case.Not, Control.Not), nrow = 2, dimnames = list(c('Case', 'Control'), c('Yes', 'No')))
      Fishers.Fit = fisher.test(Mat.Carrier)
      #allele=gsub('X', '',allele)
      outPut=data.frame(row.names = NULL, HLA=hla, allele, 
                      rsid= allele,
                      F.Pval = Fishers.Fit$p.value, F.OR = Fishers.Fit$estimate, F.CI=paste0( signif(Fishers.Fit$conf.int[1], 3),'-', signif(Fishers.Fit$conf.int[2],3)), 
                      glm.P=fitDF$p.value, glm.E=fitDF$estimate, glm.StdE = fitDF$std.error,
                      Case.Car, Case.Freq=signif(Case.Car/n.cases,3), Control.Car, 
                      Control.Freq=signif(Control.Car/n.controls, 3),
                      Case.Not, Control.Not, Case.Alleles, Case.Alelles.Not, Control.Alleles, Control.Alelles.Not, 
                      Case.Freq.Allele = signif(Case.Alleles/(2*n.cases),3), 
                      Control.Freq.Allele = signif(Control.Alleles/(2*n.controls),3), totalFreq=totalCarrier/(n.cases+n.controls), totalN=(n.cases+n.controls))
      return(outPut)
    }) 
  }
}
