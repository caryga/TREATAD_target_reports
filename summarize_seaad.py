
# libs
import scanpy as sc
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import time

# summary function -- tmp = sc.get.obs_df(fdata, [i, g, splitby])
def summarize(g, tmp, md):   
    tmp["fraction_expressed"] = tmp[g] > 0
    fraction = tmp.loc[:, ["fraction_expressed", md, cell]].groupby([md, cell]).mean()
    expression = tmp.loc[tmp[g] != 0, [g, md, cell]].groupby([md, cell]).mean().fillna(0) # mean_expression
    tmp = pd.concat([fraction, expression], axis=1).reset_index()
    tmp = tmp.rename(columns = {g:'mean_exp'})
    tmp['gene']=g
    return tmp

# read data, filter out cells from young, non-demented reference samples
st = time.time()
adata=sc.read_h5ad('SEAAD_MTG_RNAseq_final-nuclei.2022-08-18.h5ad')
adata = adata[adata.obs['Cognitive Status']!='Reference']
print('loaded adata in ',time.time()-st)
print('filtered data: ',adata.shape)

# specify dimensions by which to summarize the data
cell = 'Subclass' # broad cell classes (e.g. Astro, Micro, Sst...)
gl = adata.var.gene_ids
cog = "Cognitive Status" 
adnc = "Overall AD neuropathological Change" 

# get data with obs keys "Supertype" and either "Cognitive Status" or "ADNC"
data_obs = sc.get.obs_df(adata, [cell, cog, adnc])
data = data_obs

# loop over genes in gl and summarize expression
st1 = time.time()
for i in range(0,len(gl),500):
    
    lst = gl[i:i+499]
    
    st = time.time()
    
    tmp = sc.get.obs_df(adata, [*lst])
    tmp = pd.concat([data_obs, tmp], axis=1).reset_index()
    
    df = Parallel(n_jobs=25)(delayed(summarize)(g,tmp,cog) for g in tmp.columns[4:tmp.shape[1]] )
    
    df = pd.concat(df)
    data = pd.concat([data, df])
    
    data.to_csv('seaad_by_cog.csv')
    
    et = time.time()
    print('batch complete: ', et-st)
    

et = time.time()
print('total processing time: ', et-st1)

data.to_csv('seaad_by_cog.csv')

# EOF
