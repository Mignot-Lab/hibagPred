## author : ambati@stanford.edu
## convert DR1 alleles to amino acid calls and run associations
arg=commandArgs(trailingOnly = T)
if (length(arg) == 3){
  ## read arguments
  hlaPattern = arg[1]
  metaFile = arg[2]
  outFile = paste0(arg[3], '.csv')
  timestamp()
  message('3 ARGUMENTS SPECIFIED AS')
  message(paste0('hlaPattern :', hlaPattern, '\n', 'metaFile :', metaFile, '\n', 'outFile :', outFile))
} else {
  stop('INVALID ARGUMENTS, EXPECTING ATLEAST 3 ARGUMENTS')
}

## load the source functions
if (file.exists('scripts/funcs.R')){
  message('LOADING SOURCE FUNCTIONS')
  source('scripts/funcs.R')
} else {
  stop('SOURCE FUNCTION FILE scripts/funcs.R IS NOT READABLE')
}

## load the libs
packList = rownames(installed.packages())
libLoad=packList[grep('tidyverse|data.table|HIBAG|broom|haplo.stats', packList)]
if(length(libLoad) == 5){
  lapply(libLoad, require, character.only = T)
  message(paste0('LIBS LOADED ', libLoad, timestamp(), collapse = '\n'))
} else {
    
  message('CHECK IF LIBRARIES HIBAG, tidyverse, broom, haplo.stats & data.table ARE INSTALLED? ', timestamp())
}


## convert to HLA amino acids from HLA -DR calls
hla2aaCombine=function(hlaPattern){
  hlaDFlist = fileIO(filePre = hlaPattern)
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
    })
    names(convList) = classII
    return(convList)
  }
}

## setup asscociations
convList = hla2aaCombine(hlaPattern)
hla.meta = metaIO(metaFile = metaFile)
hlaNames = names(convList)


## run the associations
out.allele.Results=lapply(hlaNames, FUN = function(hla) 
  allele.Glm.Fishers(
  hla.meta = hla.meta, HLA_count=convList[[hla]],hla=hla)
  %>% rbindlist()) %>% rbindlist() %>% arrange(glm.P)

## write out the files
if (!is.null(out.allele.Results)){
  fwrite(out.allele.Results, file = paste0('hlaAssocs/', outFile))
  message(paste0("ASSOCIATIONS WRITTEN TO ", paste0('hlaAssocs/', outFile), " ", timestamp()))
} else {
  stop('ASSOCIATIONS FAILED')
}
#sink()








