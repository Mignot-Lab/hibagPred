## author : ambati@stanford.edu
## V3 added a single function to handle prediction and output writing, removed support for parallel preds
arg=commandArgs(trailingOnly = T)

#arg1 = model
#arg2 = plink bed files in hg19 coordinates
#arg3 = output file

#'''load the libraries or install them'''
sink(file= paste0(arg[3], '.log'))
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


## load the models

model.list <- get(load(arg[1])) 
yourgeno = hlaBED2Geno(bed.fn = paste0(arg[2], ".bed"),
                       bim.fn = paste0(arg[2], ".bim"),
                       fam.fn = paste0(arg[2], ".fam"), rm.invalid.allele = T)
summary(yourgeno)
output_pred_QC=list()
hlas_id = names(model.list)

#'''predict function to wrap in a lapply call'''
pred_hibag=function(hla, yourgeno, model.list, filePre){
  model = hlaModelFromObj(model.list[[hla]])
  print(model)
  #pred.guess = predict(model, yourgeno, type = "response+prob", cl=cl, verbose = T, match.type="Position")
  pred.guess = predict(model, yourgeno, type = "response+prob", verbose = T, match.type="Position")
  cat('PREDICTED ', hla, timestamp())
  hla.names = c(paste0(hla, '.1'), paste0(hla, '.2'), paste0(hla, '.prob'))
  temp = pred.guess$value
  names(temp)[2:4] = hla.names
  fileOut = paste0('hlaOut/', filePre,'_IMP_', hla,'.csv')
  cat('WRITING TO FILE', fileOut, timestamp())
  fwrite(temp, file=fileOut)
  #return(pred.guess$value)
}
## predict and write out the predictions

output_pred_QC=lapply(hlas_id, pred_hibag,  yourgeno = yourgeno, model.list = model.list, filePre=arg[3])
sink()

