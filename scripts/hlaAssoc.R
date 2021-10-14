## author : ambati@stanford.edu
## V1 
arg=commandArgs(trailingOnly = T)
# arg[1] = "1000g_MHC_Preds_IMP"
# arg[2] = "meta/hlaMeta.csv"
# arg[3] = "easHaps"
hlaPattern = arg[1]
metaFile = arg[2]
outFile = paste0(arg[3], '.csv')
#sink(file= paste0(arg[3], '.log'))

## load libraries
packList = rownames(installed.packages())
libLoad=packList[grep('tidyverse|data.table|HIBAG|broom|haplo.stats', packList)]
if(length(libLoad) == 5){
  lapply(libLoad, require, character.only = T)
  message(paste0('LIBS LOADED ', libLoad, timestamp(), collapse = '\n'))
} else {
    
  message('CHECK IF LIBRARIES HIBAG, tidyverse, broom, haplo.stats & data.table ARE INSTALLED? ', timestamp())
}

## functions to data IN and processing
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
    ## plot the PCS
    #pcPlot = ggplot(metaIn, aes(PC1, PC2, color=factor(Pheno)))+geom_point()
    #plotName = paste0('hlaAssocs/',arg[3], '.png')
    #message(paste0('PLOTTING PCS AND SAVING TO ', plotName, ' ', timestamp()))
    #ggsave(filename = plotName, device = 'png', height = 6.7, width = 7.7, dpi = 400)
    return(metaIn)
  } else {
    stop('META FILE DOESNT CONTAIN PHENO NAME INPUT')
  }

}

## filter hibag and one hot encode the HLA object
hibagFilter=function(hla.allele, probCutoff){
  filtered_HLA_gene = hla.allele$value[hla.allele$value$prob >= probCutoff,]
  #meanProb = mean(filtered_HLA_gene$prob)
  HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
  HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
  names(HLA_allele1) <- c("ID","allele")
  names(HLA_allele2) <- c("ID","allele")
  HLA_allele <- rbind(HLA_allele1, HLA_allele2)
  HLA_count <- as.data.frame.matrix(table(HLA_allele))
  HLA_allele_name <- make.names(names(HLA_count))
  names(HLA_count) = HLA_allele_name
  return(HLA_count)
}


#### setup associations 
allele.Glm.Fishers = function(hla, hla.meta, hla.allele, probCutoff){
  HLA_count = hibagFilter(hla.allele, probCutoff)
  HLA_allele_name = names(HLA_count)
  if (length(HLA_allele_name) >= 1){ ## more than one allele atleast
    setDT(HLA_count, key='rn', keep.rownames = T)
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
      allele=gsub('X', '',allele)
      outPut=data.frame(row.names = NULL, HLA=hla, allele, 
                      rsid= paste0(hla, ':', allele),
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


## call the functions and load the data
hlaDFlist = fileIO(filePre = hlaPattern)
hlaClassIIHaps = makeHaps(hlaDFlist = hlaDFlist)
if (exists("hlaClassIIHaps")){
  hlaDFlist$hap = hlaClassIIHaps$haplCalls
  message('HAPLOTYPE ESTIMATION SUCCESSFUL ADDING HAPLOTYPES FOR ASSOCIATIONS ', timestamp())
} else {
  message('HAPLOTYPE ESTIMATION FAILED', timestamp())
}
hlaNames = names(hlaDFlist)
hla.meta = metaIO(metaFile = metaFile)

##### make hibag object
hlaNames = names(hlaDFlist)
hla.object.list=lapply(hlaNames, function(hla){
  hla.pred.df = hlaDFlist[[hla]]
  hla.object= hlaAllele(hla.pred.df$sample.id, H1=hla.pred.df[[paste0(hla,'.1')]], 
                        H2=hla.pred.df[[paste0(hla,'.2')]], locus=hla, assembly="hg19", 
                        prob=hla.pred.df[[paste0(hla,'.prob')]])
  return(hla.object)
})
names(hla.object.list) = hlaNames

## run the associations
out.allele.Results=lapply(hlaNames, FUN = function(hla) 
  allele.Glm.Fishers(
  hla.meta = hla.meta, hla.allele=hla.object.list[[hla]],hla=hla, probCutoff = 0.3)
  %>% rbindlist()) %>% rbindlist() %>% arrange(glm.P)
## write out the files
fwrite(out.allele.Results, file = paste0('hlaAssocs/', outFile))
message(paste0("ASSOCIATIONS WRITTEN TO ", paste0('hlaAssocs/', outFile), " ", timestamp()))
#sink()

