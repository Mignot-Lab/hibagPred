from email import header
import numpy as np 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 
import os

def main():
    
    # Opts
    path = '/labs/mignot/raw_data/gwas/stanford/PMRA_PLATES_160_163/hibagPred/hlaOut/'
    prefix = 'PMRA_PLATES_160_to_163'
    
    # Parse files 
    list_files = [f for f in os.listdir(path) if prefix in f]
    list_files.sort()
    
    # First iteration 
    hla_df = pd.read_csv(os.path.join(path, list_files[0]))
    probs_df = hla_df.iloc[:,[0,3]]
    hla_df = hla_df.iloc[:,0:3]
    
    # Reamining iterations 
    for f in list_files[1:]:
        f_df = pd.read_csv(os.path.join(path, f))
        hla_df = pd.merge(left=hla_df, right=f_df.iloc[:,0:3], on='sample.id')
        probs_df = pd.merge(left=probs_df, right=f_df.iloc[:,[0,3]], on='sample.id')
        
    # Plot probs 
    plot_df = probs_df.drop('sample.id', axis=1)
    plot_df.plot.hist(subplots=True, layout=(3,4), figsize=(20,15))
    plt.show()
    plt.savefig(os.path.join(path,'PMRA_PLATES_160_163_probs_hist.png'))
    
    # Write 
    hla_df.to_csv(os.path.join(path, prefix +'_HLA.csv'), header=True, index=False)
    probs_df.to_csv(os.path.join(path, prefix + '_probs.csv'), header=True, index=False)

if __name__ == "__main__":
    main()