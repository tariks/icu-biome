# %%
import numpy as np
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins import feature_table as ft
from qiime2.plugins import taxa
# %%
asv = Artifact.load('./icu-biome/qiime/table_stool.qza')
tax = Artifact.load('./icu-biome/qiime/tax_curated.qza')
# %%
#filt_asv = pd.read_csv('./icu-biome/feature_tables/59N_ASV_rare6k_filt002p7.csv',index_col=0).columns.to_list()
meta = pd.read_csv('./icu-biome/meta/59N_alpha.csv',index_col=0)
# %%
ftable = ft.actions.filter_features_conditionally(asv,.00002,.01).filtered_table
ftable = ft.actions.filter_samples(ftable,5000).filtered_table
# %%
biggenus = taxa.actions.collapse(ftable,tax,7).collapsed_table.view(pd.DataFrame)
biggenus.to_csv('./icu-biome/feature_tables/59N_genus.csv')
ftable = taxa.actions.collapse(ftable,tax,7).collapsed_table
biggenus.shape

# %%
genus = ft.actions.filter_features_conditionally(ftable,0.000000000001,.2).filtered_table.view(pd.DataFrame)
genus = genus.loc[genus.sum(axis=1)>1000]

# %%
def parseTaxonomy(table,thres=60):
    # extracts genus and applies IDTaxa threshold
    # returns clean matrix
    df = pd.DataFrame(index=table.index)
    cols=[]
    for i in table.columns:
        if 'genus' in i or 'Enterobacteriaceae' in i:
            cols.append(i)
    x=table[cols].copy()
    for i in cols:
        s = i.split(');')
        if '__' in i:
            s=s[-2]
        else:
            s=s[-1]
        name = s.split(' (')[0]
        s = s.split(' ')
        if int(s[-1].split('%')[0]) > thres:
            if name in df.columns:
                df[name]+=table[i]
            else:
                df[name]=table[i].copy()
    return df.loc[:,df.sum()>0].astype(int)

# %%
x=parseTaxonomy(genus,thres=70)
x.to_csv('./icu-biome/feature_tables/59N_genus_t70.csv')
# %%
x.loc[meta['Day']==1].to_csv('./icu-biome/feature_tables/52N_genus_t70.csv')

# %%
from skbio.stats.composition import multiplicative_replacement
# %%
x[:]=multiplicative_replacement(x)
x.loc[meta['Day']==1].to_csv('./icu-biome/feature_tables/52N_genus_t70_nozero.csv')
x[:]=np.log(x)
x.to_csv('./icu-biome/feature_tables/59N_genus_t70_log.csv')
# %%
meta.loc[meta['Day']==1,['APACHE','Age']].to_csv('./icu-biome/selbal/covariates.csv')
meta.loc[meta['Day']==1,'month'].astype(int).to_csv('./icu-biome/selbal/target.csv')
# %%
