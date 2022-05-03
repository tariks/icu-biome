# %%
from qiime2 import Artifact
import numpy as np
import pandas as pd
import seaborn as sns
from skbio import diversity as div

# %%
#T = pd.read_csv("./icu-biome/feature_tables/all_asv.csv", index_col=0)
T = Artifact.load('./icu-biome/qiime/table_stool.qza').view(pd.DataFrame)
meta = pd.read_csv("./icu-biome/meta/59meta.csv", index_col=0)

T = T.loc[T.index.isin(meta.index)]
#T.loc[meta['Day']==1].T.sum().describe()

# %%
def rarefaction(M, seed=0, depth=6000):
    prng = np.random.RandomState(seed)  # reproducible results
    noccur = np.sum(M, axis=1)  # number of occurrences for each sample
    nvar = M.shape[1]  # number of variables

    Mrarefied = np.empty_like(M)
    for i in range(M.shape[0]):  # for each sample
        p = M[i] / float(noccur[i])  # relative frequency / probability
        choice = prng.choice(nvar, depth, p=p)
        Mrarefied[i] = np.bincount(choice, minlength=nvar)

    return Mrarefied


# %%
T[:] = rarefaction(T.values)
# %%
T = T.loc[:,T.sum()>7] # .002% threshold. 8k * 51 * .00002 = 8.16
T.sum().describe()
# %%
div_metrics= [
    "ace",
    "chao1",
    "dominance",
    "fisher_alpha",
    "pielou_e",
    "shannon",
    "simpson",
]
for i in div_metrics:
    meta[i] = div.alpha_diversity(i, T.astype(int).values, ids=T.index,)
    print(meta[i].describe().T)
# %%
%matplotlib inline
sns.pairplot(data=meta,vars=['Death']+div_metrics,hue='Enterotype')

# %%
#meta.to_csv('./icu-biome/meta/59N_alpha.csv')
meta.loc[meta['Day']==1].to_csv('./icu-biome/meta/52_alpha.csv')
# %%
T.to_csv('./icu-biome/feature_tables/59N_ASV_rare6k_filt002p7.csv')
# %%
