# %%
%matplotlib inline
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
gpac = [
    'Anaerococcus',
    'Fenollaria',
    'Finegoldia',
    'Peptococcus',
    'Peptostreptococcus',
    'Coprococcus',
    'Atopobium',
    'Ruminococcus',
    'Parvimonas',
    'Peptoniphilus',
    'Blautia',
    'Gallicola',
    'Murdochiella',
    'Sarcina',
]
mypal= [
(0.577, 0.556, 0.324, 1.000),
(0.749, 0.647, 0.440, 1.000),
(0.941, 0.808, 0.643, 1.000),
(0.998, 0.963, 0.797, 1.000),

(0.129, 0.314, 0.467, 1.000),
(0.285, 0.555, 0.761, 1.000),
(0.438, 0.661, 0.724, 1.000),
(0.623, 0.808, 0.866, 1.000),

(0.675, 0.396, 0.569, 1.000),
(0.816, 0.527, 0.669, 1.000),
(0.831, 0.628, 0.739, 1.000),
(0.889, 0.721, 0.839, 1.000),

]
mypal = sns.color_palette(mypal,desat=.9)
# %%
plt.style.use("../lance.txt")
meta = pd.read_csv("../meta/52_bal.csv", index_col=0)

phygen = {
    'Bacteroidetes': ['Bacteroides','Parabacteroides','Alistipes','Other Bacteroidetes'],
    'Proteobacteria': ['Enterobacteriaceae','Pseudomonas','Acinetobacter','Other Proteobacteria'],
    'Firmicutes': ['GPAC','Enterococcus','Staphylococcus','Other Firmicutes'],
}
colors = {
#    'Firmicutes': sns.color_palette('tab20b',4)[-4:][::-1],
#    'Proteobacteria': sns.color_palette('tab20b',12)[-4:][::-1], 
#    'Bacteroidetes': sns.color_palette('tab20c',4),
    'Firmicutes': mypal[:4][::-1],
    'Proteobacteria': mypal[-4:][::-1],
    'Bacteroidetes': mypal[4:8][::-1],
}
taxa = []
c = []
for i in ['Firmicutes','Proteobacteria','Bacteroidetes']:
    taxa+=phygen.get(i)[::-1]
    c+=list(colors.get(i))
# %%
phylum = pd.read_csv("../feature_tables/52_phyla70.csv", index_col=0)
genus = pd.read_csv("../feature_tables/52N_genus_t70.csv", index_col=0)
genus['GPAC']=0
for i in genus.columns:
    if i in gpac:
        genus['GPAC']+=genus[i]
        del genus[i]
phylum['Bacteroidetes'] = phylum['Bacteroidota'].copy()
phylum = phylum[['Bacteroidetes','Proteobacteria','Firmicutes']]
for i in ['Bacteroidetes','Proteobacteria','Firmicutes']: genus['Other '+i]=0
#genus = genus.divide(genus.sum(axis=1),axis=0)
genus['Other Firmicutes'] = phylum['Firmicutes']-genus[phygen['Firmicutes']].sum(axis=1)
genus['Other Proteobacteria'] = phylum['Proteobacteria']-genus[phygen['Proteobacteria']].sum(axis=1)
genus['Other Bacteroidetes'] = phylum['Bacteroidetes']-genus[phygen['Bacteroidetes']].sum(axis=1)
print(genus.mean().sort_values(ascending=False).index[:20].to_list())

genus = genus[taxa]
#phylum = phylum.divide(phylum.sum(axis=1),axis=0)
phylum[:] = multiplicative_replacement(phylum)
genus[:] = multiplicative_replacement(genus)
#genus = genus.divide(genus.sum(axis=1),axis=0)

order = meta['mbal'].sort_values(ascending=True).index
genus=genus.loc[order]

params=dict(
        width=1,
        lw=.25,
        ec='#ababab',
        snap=True,
#        alpha=.9,
)
# %%
h= np.ones(52)
fig,ax = plt.subplots()
fig.set(figheight=2.5,figwidth=4)
genus=genus.loc[phylum['Bacteroidetes'].sort_values()[::-1].index]
genus=genus.loc[genus['Bacteroides'].sort_values()[::-1].index]
for j,i in enumerate(taxa[::-1]):
    ax.bar(x= np.arange(1,53),
        height=h,
#        bottom=h-genus[i].values,
        color=c[-1-j],
        label=i,
        **params,
    )
    h-=genus[i].values
ax.set_xlabel('')
ax.set_xticklabels([])
ax.set_ylim(0,1)
ax.set_xlim(.5,52.5)
ax.grid(False)
ax.yaxis.set_major_formatter(ticker.PercentFormatter(1))
ax.legend(loc='upper left',bbox_to_anchor=[1,1],borderpad=0,)
#ax.set_xticks([1,13.5,26.5,39.5,52])
#ax.set_xticklabels([-4.07,.57,1.9,3.97,6.71])
ax.set_title('Composition of stool bacteria in 52 patients at time of ICU admission',loc='center')
plt.savefig('../plots/taxabars.pdf',dpi=300,transparent=True,bbox_inches='tight')


# %%
