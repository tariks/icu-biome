# %%
from qiime2 import Artifact
import numpy as np
import pandas as pd
import seaborn as sns
from skbio import diversity as div

# %%
# table = Artifact.load('../qiime/table_stool.qza').view(pd.DataFrame)
table = Artifact.load("../qiime/feature-table2.qza").view(pd.DataFrame)
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)

table = table.loc[table.index.isin(meta.index)]
table

# %%
def rarefaction(M, seed=0, depth=6700):
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
conv = {
    "ace": "ACE (R)",
    "chao1": "Chao1 (R)",
    "dominance": "Dominance",
    "fisher_alpha": "Fisher (D)",
    "pielou_e": "Pielou (E)",
    "shannon": "Shannon (D)",
}
# %%
# T = T[[i for i in T.columns if (T[i]>1).sum()>1]]
T = table.copy()
T[:] = rarefaction(table.values, depth=8000)
# rel = T.loc[:,T.sum()>7] # .002% threshold. 8k * 51 * .00002 = 8.16
# T = T[[i for i in T.columns if (T[i]>0).sum()>1]]
T = T[[i for i in T.columns if T[i].sum() > 20]]
T.sum(axis=1).describe()
# %%
div_metrics = [
    "ace",
    "chao1",
    "fisher_alpha",
    "pielou_e",
    "shannon",
    "dominance",
]
for i in div_metrics:
    meta[conv.get(i)] = div.alpha_diversity(
        i,
        T.astype(int).values,
        ids=T.index,
    )
    print(meta[conv.get(i)].describe().T)
# %%
meta["Simpson (D)"] = 1 / meta["Dominance"]
meta["Simpson (E)"] = div.alpha_diversity(
    "simpson_e", T.astype(int).values, ids=T.index
)

# %%
# meta.to_csv('./icu-biome/meta/59N_alpha.csv')
meta.to_csv("../meta/meta52_current.csv")
# %%
T.to_csv("./icu-biome/feature_tables/59N_ASV_rare6k_filt002p7.csv")
# %%
