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
(0.129, 0.314, 0.467, 1.000),
(0.200, 0.518, 0.761, 1.000),
(0.400, 0.706, 0.792, 1.000),
(0.610, 0.804, 0.862, 1.000),
(0.337, 0.318, 0.200, 1.000),
(0.553, 0.522, 0.337, 1.000),
(0.784, 0.761, 0.565, 1.000),
(0.911, 0.896, 0.756, 1.000),
(0.995, 0.995, 0.819, 1.000),
(0.984, 0.922, 0.600, 1.000),
(0.953, 0.800, 0.404, 1.000),
(0.922, 0.655, 0.329, 1.000),

]
mypal = sns.color_palette(mypal,desat=1)
# %%
plt.style.use("../lance.txt")
meta = pd.read_csv("../meta/52_bal.csv", index_col=0)

phygen = {
    'Bacteroidetes': ['Bacteroides','Parabacteroides','Porphyromonas','Other Bacteroidetes'],
    'Proteobacteria': ['Enterobacteriaceae','Pseudomonas','Acinetobacter','Other Proteobacteria'],
    'Firmicutes': ['GPAC','Enterococcus','Staphylococcus','Other Firmicutes'],
}
colors = {
#    'Firmicutes': sns.color_palette('tab20b',4)[-4:][::-1],
#    'Proteobacteria': sns.color_palette('tab20b',12)[-4:][::-1], 
#    'Bacteroidetes': sns.color_palette('tab20c',4),
    'Firmicutes': mypal[4:8][::-1],
    'Proteobacteria': mypal[-4:],
    'Bacteroidetes': mypal[:4][::-1],
}
taxa = []
c = []
for i in ['Proteobacteria','Firmicutes','Bacteroidetes']:
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
#genus = genus.divide(genus.sum(axis=1),axis=0)
print(genus.max().sort_values(ascending=False).index[:20].to_list())

genus = genus[taxa]
#phylum = phylum.divide(phylum.sum(axis=1),axis=0)
phylum[:] = multiplicative_replacement(phylum)
genus[:] = multiplicative_replacement(genus)

order = meta['mbal'].sort_values(ascending=True).index
genus=genus.loc[order]

params=dict(
        width=1,
        lw=.2,
        ec='#919191',
        snap=True,
#        alpha=.9,
)
# %%
h= np.zeros(52)
fig,ax = plt.subplots()
fig.set(figheight=2.5,figwidth=4)
phylum=phylum.sort_values(by=['Bacteroidetes','Firmicutes'],ascending=False)
genus=genus.loc[phylum.index]
for j,i in enumerate(taxa[::-1]):
    ax.bar(x= np.arange(1,53),
        #height=h,
        height=genus[i],
        #bottom=h-genus[i].values,
        bottom=h,
        color=c[-1-j],
        label=i,
        **params,
    )
    h+=genus[i].values
ax.set_xlabel('')
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_ylim(0,1)
ax.set_xlim(.5,52.5)
ax.grid(False)
#ax.set_xticks([1,13.5,26.5,39.5,52])
#ax.set_xticklabels([-4.07,.57,1.9,3.97,6.71])
ax.set_title('Composition of stool bacteria in 52 patients at time of ICU admission',loc='center')
ax.invert_yaxis()
h,l=ax.get_legend_handles_labels()
h,l=h[::-1],l[::-1]
ax.legend(loc='upper left',bbox_to_anchor=[1,1],borderpad=0,)
ax.set_yticklabels(['0%','20%','40%','60%','80%','100%'][::-1])
#ax.yaxis.set_major_formatter(ticker.PercentFormatter(1))


plt.savefig('../plots/taxabars.pdf',dpi=300,transparent=True,bbox_inches='tight')


# %%
