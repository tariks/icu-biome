# %%
import numpy as np
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins import feature_table as ft
from qiime2.plugins import taxa

# %%
asv = Artifact.load("../qiime/table_stool.qza")
tax = Artifact.load("../qiime/tax_curated.qza")
# %%
meta = pd.read_csv("../meta/52_alpha.csv", index_col=0)
# %%
# ftable = ft.actions.filter_features_conditionally(asv,.001,.01).filtered_table
ftable = ft.actions.filter_samples(asv, 1000).filtered_table
# %%
ftable = taxa.actions.collapse(ftable, tax, 7).collapsed_table
ftable = ftable.view(pd.DataFrame)
ftable.shape

# %%
def freq(table):
    # return taxa proportion of total
    return table.sum() / table.sum().sum()


def freq2(table):
    # return prevalence across cohort
    return (table > 0).astype(int).sum() / table.shape[0]


def freqprint(table):
    s1, s2 = freq(table), freq2(table)
    s2 = s2.sort_values()[::-1]
    for i in s1.index:
        print("{}\t{:.2f}\t{:.2f}".format(i, s1[i], s2[i]))
    return s2


# %%
import re
from collections import defaultdict

phyladict = defaultdict(str)


def parseTaxonomy(table, thres=60):
    # extracts genus and applies IDTaxa threshold
    df = pd.DataFrame(index=table.index)
    cols = table.columns.to_list()
    cols = [i for i in cols if "%" in i]
    cols = [i for i in cols if "phylum" in i]
    #    cols = [i for i in cols if "genus" in i or "Enterobacteriaceae" in i]
    for i in cols:
        s = i.split(";")
        phyla = s[2].split()[0]
        for j in s[::-1]:
            if "Enterobacteriaceae" in j:
                s = j
                break
            elif "genus" in j:
                s = j
                break
            elif "phylum" in j:
                s = j
                break
            else:
                pass
        name = re.split("\s\(", s)[0]
        phyladict[name] = phyla
        score = int(re.findall("\d+%", s)[0][:-1])
        if score > thres:
            if name in df.columns:
                df[name] += table[i]
            else:
                df[name] = table[i].copy()
    df = df.loc[:, (df.sum() > 0)]
    return df.astype(int)


# %%
big = parseTaxonomy(ftable, thres=70)
big = big.loc[meta.index]
s2 = freqprint(big).sort_values(ascending=False)
print(big.shape)
# %%
p = pd.Series(phyladict)
p[p == "Bacteroidota"] = "Bacteroidetes"
p.to_csv("../feature_tables/taxmap.csv", header=False)

# %%
# preview matrix
x = big.loc[meta.index]
print("Prevalence:")
for i in [0.2, 0.25, 0.3]:
    print(">{}\t{}".format(str(i)[1:], (s2 > i).sum()))
print("Samples < 5:")
for i in [0.2, 0.25, 0.3]:
    # cnt = (x.loc[:,s2>i].sum(axis=1) < 1000).sum()
    cnt = ((x.loc[:, s2 > i] > 0).sum(axis=1) < 10).sum()
    print(">{}\t{}".format(str(i)[1:], cnt))


# %%
from skbio.stats.composition import multiplicative_replacement

# %%
idx = s2 > 0.2
idxrow = (x.loc[:, idx] > 0).sum(axis=1) >= 0
x.loc[idxrow, idx].to_csv("../feature_tables/52N_genus_t70.csv")
x[:] = multiplicative_replacement(x)
x.loc[idxrow, idx].to_csv("../feature_tables/52N_genus_t70_nozero.csv")
meta.loc[idxrow, "month"].to_csv("../selbal/target.csv")
meta.loc[idxrow, ["APACHE", "Age"]].to_csv("../selbal/covariates.csv")
# %%
