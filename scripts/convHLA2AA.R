## author : aditya.ambati@gmail.com
## convert class2 alleles to amino acid calls dosages including haplotypes for DRB1
arg=commandArgs(trailingOnly = T)
if (length(arg) == 4){
  ## read arguments
  DQB1 = arg[1]
  DQA1 = arg[2]
  DRB1 = arg[3]
  outFile = arg[4]
  fileList = c(arg[1:3])
  timestamp()
  message('4 ARGUMENTS SPECIFIED AS')
  message(paste0('DQB1 :', DQB1, '\n', 'DQA1 :', DQA1, '\n', 'DRB1 :', DRB1, '\n','outFile :', outFile))
} else {
  stop('INVALID ARGUMENTS, EXPECTING ATLEAST 3 ARGUMENTS')
}

packList = rownames(installed.packages())
libLoad=packList[grep('tidyverse|data.table|HIBAG|broom|haplo.stats|stringr', packList)]
if(length(libLoad) == 6){
  lapply(libLoad, require, character.only = T)
  message(paste0('LIBS LOADED ', libLoad, timestamp(), collapse = '\n'))
} else {
  
  message('CHECK IF LIBRARIES HIBAG, tidyverse, broom, haplo.stats & data.table ARE INSTALLED? ', timestamp())
}

# load the source functions
if (file.exists('scripts/funcs.R')){
  message('LOADING SOURCE FUNCTIONS')
  #source('scripts/HLA2AA.R')
  source('scripts/funcs.R')
} else {
  stop('SOURCE FUNCTION FILE  scripts/funcs.R IS NOT READABLE')
}

if (!dir.exists('aminoAcids/')){
  dir.create('aminoAcids/')
}

fileIO=function(fileList){
  pathTohla=c(fileList)#list.files(path = 'hlaOut/', pattern = paste0(filePre, '*'), full.names = T)
  if (length(pathTohla) >= 3){
    message(paste0('IDENTIFIED HLA IMPUTE FILES IN hlaOut/ ', pathTohla, ' ', timestamp(), collapse = '\n'), appendLF = T)
    hlaNames = c('DQB1', 'DQA1', 'DRB1')#gsub('IMP_|.csv', '', unique(str_extract(pattern = "DQA1|DRB1|DQB1", pathTohla)))
    print(pathTohla)
    hlaDFlist = lapply(pathTohla, function(hlaFile){
      temp = fread(hlaFile, key='sample.id', select = c(1:4)) # take first 4 cols to output
    })
    names(hlaDFlist) = hlaNames
    return(hlaDFlist)
  } else {
    message('HAVE IMPUTATIONS BEEN PERFORMED ? RUN HibagPred.R FIRST OR CHECK DIR LOCATION ', timestamp() )
    stop()
  }
}

hla2aaCombine=function(fileList, outFile){
  hlaDFlist = fileIO(fileList)
  classII = c('DQB1', 'DQA1', 'DRB1')
  if(!any(names(hlaDFlist) %in% classII)){
    stop('INVALID LOCI SPECIFIED PLEASE CHECK IF DQB1 DQA1 DRB1 EXIST in hlaOut FOLDER')
  } else {
    convList=lapply(classII, function(hla){
      #make sub arguments to pass HLA2AA function
      hlaFile = hlaDFlist[[hla]]
      locus = hla
      out=HLA2AA(hlaFile,locus, probCutoff=0.3)
      convDF = data.table(sample.id=rownames(out[[1]]))
      for(pos in out){
        temp=data.table(pos, keep.rownames = T)
        if(identical(temp$rn, convDF$sample.id)){
          convDF = cbind.data.frame(convDF, pos)
        }
      }
      if(hla == 'DRB1'){
        drb1aalist=aminoDRB1(hlaFile = hlaFile, locus = 'DRB1', probCutoff = 0.3)
        haps=aminoHap(drb1aalist)
        print(haps)
        print(names(convDF))
        if(identical(haps$rn, convDF$sample.id)){
          names(haps) = paste0(hla,'_', names(haps))
          convDF = cbind.data.frame(convDF, haps[, -1])
        } else {
          stop()
        }
        
      }
      setDT(convDF)
      fwrite(convDF, file = paste0('aminoAcids/', outFile,'_', hla, '_AA.csv'))
      print(dim(convDF))
    })
    
    names(convList) = classII
    return(convList)
  }
}


#convList = hla2aaCombine(hlaPattern="PPMI_Euro")
#fileList = c('hlaOut//1000g_MHC_Preds_IMP_DQB1.csv', 'hlaOut//1000g_MHC_Preds_IMP_DQA1.csv', 'hlaOut//1000g_MHC_Preds_IMP_DRB1.csv')#apply(str_split(list.files('hlaOut/', pattern = '1000g_MHC_Preds_IMP_*'), pattern = "_", simplify = T)[, 1:4], 1, paste0, collapse = "_") %>% unique()
print('RUNNING THE PIPELINE')
timestamp()
hla2aaCombine(fileList, outFile)
#e.g run 
#Rscript scripts/convHLA2AA.R hlaOut//1000g_MHC_Preds_IMP_DQB1.csv hlaOut/1000g_MHC_Preds_IMP_DQA1.csv hlaOut/1000g_MHC_Preds_IMP_DRB1.csv 1000gtest
