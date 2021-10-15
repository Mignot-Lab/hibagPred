# HibagPred.R
A script to predict “A”    “B”    “C”    “DPA1" “DPB1” “DQA1" “DQB1” “DRB1" “DRB3” “DRB4" “DRB5” at 4-digit resolution using HIBAG.

Best to use imputed HLA region (preferably 1000 genomes phase3 or phase1, chr6:29519561-33083508) as input in plink binary format to the script. Make sure there are no duplicate positions or rsids before you pass the plink files to the script.

A typical usage would require 3 arguments (below).

Argument 1 = the trained model in this case (ask ambati@stanford.edu for models).

Argument 2 = is the plink binary file that contains genotypes (preferably imputed ones to increase coverage with the snps in model) (prefix only).

Argument 3 = is the path to the output filename (prefix only).

 

Example usage in linux or mac os terminal would be from your working directory containing the plink binaries and the model along with attached script.

```Rscript scripts/HibagPred.R models/model.Rdata plink_binary_prefix output_prefix```

This requires HIBAG package and data table package which the script will automatically install if it doesn’t detect them.
After the script runs, you will have the respective output file with HLA calls and probabilities  and an associated log file.

# hlaAssoc.R
A script to run associations from the predicted HLA outs in hlaOut folder
This script requires the following dependencies
1. data.table
1. HIBAG
1. tidyverse
1. broom

The script will fit a logistic regression carrier model and also compute fisher's test on carrier counts between cases and controls
The script requires that HibagPred.R is run first, the script requires 3 arguments

Argument 1 = the file prefix of the imputed HLA, in this example it is 1000g_MHC_Preds_IMP.

Argument 2 = is the meta file that must contain the  sample.id, phenotype and the first 4 genetic prinicipal components. The column names must be (sample.id,PC1,PC2,PC3,PC4,Pheno) , see example in meta/hlaMeta.csv, 
the pheno must all be binary 1 or 0 and missing can be labelled as NA

Argument 3 = is the path to the output filename (prefix only).

Example usage in linux or mac os terminal would be from your working directory.

```Rscript scripts/hlaAssoc.R EAS1KG_IMP meta/hlaMeta.csv easHap```

If the script ran successfully, you will have the associations in the folder hlaAssocs/*.csv , see example in hlaAssocs/easTest.csv

# HLA2AA.R
A script to make Amino Acid calls and run associations from HLA-Class II  (DQB1, DQA1, DRB1) predicted HLA outs in hlaOut folder
This script requires the following dependencies
1. data.table
1. HIBAG

The script requires the following arguments

Argument 1 = the file prefix of the imputed HLA, in this example it is 1000g_MHC_Preds_IMP.

Argument 2 = is the meta file that must contain the  sample.id, phenotype and the first 4 genetic prinicipal components. The column names must be (sample.id,PC1,PC2,PC3,PC4,Pheno) , see example in meta/hlaMeta.csv, 
the pheno must all be binary 1 or 0 and missing can be labelled as NA

Argument 3 = is the path to the output filename (prefix only).

Example usage in linux or mac os terminal would be from your working directory.

```Rscript scripts/HLA2AA.R EAS1KG_IMP meta/hlaMeta.csv easAA```