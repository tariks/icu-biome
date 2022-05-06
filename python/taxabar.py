# %%
#%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
from scipy import stats
from skbio.stats.composition import multiplicative_replacement, clr, ilr
from sklearn.decomposition import PCA
import cmasher as cmr

# %%
plt.style.use("../lance.txt")
meta = pd.read_csv("../meta/52_bal.csv", index_col=0)

phygen = {
    'Bacteroidetes': ['Bacteroides','Parabacteroides','Alistipes','Other Bacteroidetes'],
    'Proteobacteria': ['Enterobacteriaceae','Acinetobacter','Pseudomonas','Other Proteobacteria'],
    'Firmicutes': ['Enterococcus','Staphylococcus','Anaerococcus','Other Firmicutes'],
}
colors = {
    'Firmicutes': sns.color_palette('tab20b',4)[-4:][::-1],
    'Proteobacteria': sns.color_palette('tab20b',12)[-4:][::-1], 
    'Bacteroidetes': sns.color_palette('tab20c',4),
}
taxa = []
c = []
for i in ['Bacteroidetes','Proteobacteria','Firmicutes',]:
    taxa+=phygen.get(i)[::-1]
    c+=list(colors.get(i))
# %%
phylum = pd.read_csv("../feature_tables/52_phyla70.csv", index_col=0)
genus = pd.read_csv("../feature_tables/52N_genus_t70.csv", index_col=0)

phylum['Bacteroidetes'] = phylum['Bacteroidota'].copy()
phylum = phylum[['Bacteroidetes','Proteobacteria','Firmicutes']]
#genus = genus.divide(genus.sum(axis=1),axis=0)
genus['Other Firmicutes'] = phylum['Firmicutes']-genus[phygen['Firmicutes'][:-1]].sum(axis=1)
genus['Other Proteobacteria'] = phylum['Proteobacteria']-genus[phygen['Proteobacteria'][:-1]].sum(axis=1)
genus['Other Bacteroidetes'] = phylum['Bacteroidetes']-genus[phygen['Bacteroidetes'][:-1]].sum(axis=1)
genus = genus[taxa]
#phylum = phylum.divide(phylum.sum(axis=1),axis=0)
#genus = genus.divide(genus.sum(axis=1),axis=0)
#phylum[:] = multiplicative_replacement(phylum)
genus[:] = multiplicative_replacement(genus)

taxa

# %%
order = meta['mbal'].sort_values(ascending=True).index
genus=genus.loc[order]

# %%
params=dict(
        width=1,
        lw=.5,
        ec='w',
        snap=True,
        alpha=1,
)
# %%
x=genus
h= np.ones(52)
fig,ax = plt.subplots()
fig.set(figheight=2.5,figwidth=4)
for j,i in enumerate(x.columns[::-1]):
    ax.bar(x=np.arange(1,53),
        height=h,
        bottom=h-x[i],
        color=c[-1-j],
        label=i,
        **params,
    )
    h-=x[i]
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_ylim(0,1)
ax.set_xlim(.5,52.5)
ax.grid(False)
ax.yaxis.set_major_formatter(ticker.PercentFormatter(1))
ax.legend(loc='upper left',bbox_to_anchor=[1,1],borderpad=0,)
ax.set_xticks([1,13.5,26.5,39.5,52])
ax.set_xticklabels([-4.07,.57,1.9,3.97,6.71])
ax.set_title('Composition of stool bacteria in 52 patients at time of ICU admission',loc='center')
plt.savefig('../plots/taxabars.pdf',dpi=300,transparent=True,bbox_inches='tight')

# %%

# %%
