# hibagPred
A script to predict “A”    “B”    “C”    “DPA1" “DPB1” “DQA1" “DQB1” “DRB1" “DRB3” “DRB4" “DRB5” at 4-digit resolution using HIBAG.

Best to use imputed HLA region (preferably 1000 genomes phase3 or phase1, chr6:29519561-33083508) as input in plink binary format to the script. Make sure there are no duplicate positions or rsids before you pass the plink files to the script.

A typical usage would require 3 arguments (below).

Argument 1 = the trained model in this case (ask ambati@stanford.edu for models).

Argument 2 = is the plink binary file that contains genotypes (preferably imputed ones to increase coverage with the snps in model) (prefix only).

Argument 3 = is the path to the output filename (prefix only).

 

Example usage in linux or mac os terminal would be from your working directory containing the plink binaries and the model along with attached script.

#Rscript HibagPred.R  model.Rdata    plink_binary_prefix  output_prefix

This requires HIBAG package and data table package which the script will automatically install if it doesn’t detect them.

 

After the script runs, you will have the respective output file with HLA calls and probabilities  and an associated log file.

