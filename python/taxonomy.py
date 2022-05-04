# %%
import numpy as np
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins import feature_table as ft
from qiime2.plugins import taxa
# %%
asv = Artifact.load('../qiime/table_stool.qza')
tax = Artifact.load('../qiime/tax_curated.qza')
# %%
meta = pd.read_csv('../meta/52_alpha.csv',index_col=0)
# %%
#ftable = ft.actions.filter_features_conditionally(asv,.00002,.01).filtered_table
ftable = ft.actions.filter_samples(asv,1000).filtered_table
# %%
biggenus = taxa.actions.collapse(ftable,tax,7).collapsed_table.view(pd.DataFrame)
biggenus.to_csv('../feature_tables/52N_genus.csv')
ftable = taxa.actions.collapse(ftable,tax,7).collapsed_table
biggenus.shape


# %%
def freq(table):
    # return taxa proportion of total 
    return table.sum() / table.sum().sum()
def freq2(table):
    # return prevalence across cohort
    return (table>0).astype(int).sum() / table.shape[0]
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
big = parseTaxonomy(biggenus,thres=70)
print(big.shape)
s = freq(big).sort_values(ascending=False)
s2 = freq2(big).sort_values(ascending=False)
for i in s.index:
    print('{}\t{:.2f}\t{:.2f}'.format(i,s[i],s2[i]))
# %%
x=parseTaxonomy(biggenus,thres=70)
# %%
x.to_csv('../feature_tables/59N_genus_t70.csv')
# %%
x=x.loc[meta.index]
print(x.shape)

# %%
from skbio.stats.composition import multiplicative_replacement
# %%
x.loc[:,s>.2].to_csv('../selbal/52N_genus_t70.csv')
x[:]=multiplicative_replacement(x)
x.loc[:,s>.2].to_csv('../selbal/52N_genus_t70_nozero.csv')
# %%
