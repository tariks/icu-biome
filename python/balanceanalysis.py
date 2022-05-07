# %%
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
from skbio.stats.composition import multiplicative_replacement, clr, ilr
from sklearn.decomposition import PCA
import cmasher as cmr
from pca import pca

# %%
meta = pd.read_csv("../meta/52_bal.csv", index_col=0)
table = pd.read_csv("../feature_tables/52N_genus_t70.csv", index_col=0)
# %%
#table = pd.read_csv("../feature_tables/52_phylum80.csv", index_col=0)
s = table.sum() / table.sum().sum()
s = s.sort_values()[::-1]
s1 = (table > 0).sum().sort_values()[::-1]
a, b = set(s.index[:3].to_list()), set(s1.index[:3].to_list())
c = list(a.union(b))
table = table.loc[:, table.columns.isin(c)]
table.shape
# %%
plt.style.use("../lance.txt")
# %%
x=table.copy()
x[:]=clr(multiplicative_replacement(x))
top=x.sum().sort_values()[::-1].index[:50].to_list()
top
# %%

mod = pca(normalize=True)
results = mod.fit_transform(x[top])
fig,ax = mod.biplot(
    y = (meta['mbal']>2).astype(int),
    n_feat=3,
    PC=[0,1],
#    SPE=True,
#    hotellingt2=True,
    figsize=(2.5,2.5),
    alpha_transparency=.3,
    visible=False,
    color_arrow='#810004',
    verbose=0,
    title='PCA analysis',
    legend=False,
)
for t in ax.texts:
    if len(t.get_text())<5:
        t.remove()
    if '(' in t.get_text():
        s = t.get_text().split('(')
        number = float(s[1].split(')')[0])
        s=s[0]+'({:.2f})'.format(number)
        a,b=t.get_position()
        plt.setp(t,text=s,color='#333333',fontsize=5,x=a*1.3,y=b*1.3)
for marker in ax.collections:
    marker.set_sizes(marker.get_sizes()/8)
ax.collections[0].set_color('#810004')
ax.collections[1].set_color('#004481')
ax.set_visible(True)
fig.set_visible(True)
ax.legend(labels=['Low balance','High balance'],
    loc='upper right',
    bbox_to_anchor=[.99,1.15],
    markerscale=.9,
    handletextpad=.1,
)
plt.savefig('../pca_genus.pdf',dpi=300,transparent=True,bbox_inches='tight')
plt.show()

# %%
sns.histplot(data=meta, x="mbal", hue="month")
# %%
cont = ["Age", "APACHE", "shannon", "simpson", "Chao1", "Death"]
sns.pairplot(data=meta, x_vars="mbal", y_vars=cont)

# %%
meta["M"] = (meta["mbal"] > 1.7).astype(int)
meta.groupby("M").describe()["Death"]
# %%
sns.scatterplot(
    data=meta,
    x="mbal",
    y="Death",
    snap=True,
)
# %%
fig, ax = plt.subplots()
meta['d']=meta['Death'].copy()
meta.loc[meta['Death']>180,'d'] = 180

sns.swarmplot(data=meta, x="M", y="Death",
    palette='dark',hue='month',
   # s=2.1,
    alpha=.9,ec='w',lw=.5)

fig.set(figheight=1.5,figwidth=2)

ax.set_xlabel('')
ax.set_xticklabels([
    '<0','<1.7','>4'
])
ax.set_ylim(0,380)
ax.set_xlim(-.5,2.5)

# %%
x = table.copy()
x[:] = multiplicative_replacement(x)
#x = pd.DataFrame(ilr(x), index=x.index)
x[:] = clr(x)
pca = PCA(whiten=True)
pca.fit(x)
fig, ax = plt.subplots()
pcs = pca.fit_transform(x)
pcs = pd.DataFrame(pcs, index=x.index)
g=sns.scatterplot(
    x=pcs[0],
    y=pcs[1],
    hue=meta["mbal"],
#    palette='dark',
#    palette='flare',
#    palette=cmr.dusk,
    palette=cmr.amber,
#    palette=cmr.bubblegum,
    s=20,
    alpha=0.9,
    ax=ax,
    snap=True,
)
ax.text(
    x=0.99,
    y=0.01,
    s="PC2 {:.2%}".format(pca.explained_variance_ratio_[1]),
    ha="right",
    va="bottom",
    transform=ax.transAxes,
    fontsize=6,
)
ax.text(
    x=0.01,
    y=0.99,
    s="PC1 {:.2%}".format(pca.explained_variance_ratio_[0]),
    ha="left",
    va="top",
    transform=ax.transAxes,
    fontsize=6,
)
ax.margins(0.12)
ax.legend().set_title("")
ax.set_xlabel("")
ax.set_ylabel("")
fig.set(figheight=2.,figwidth=3)
sns.move_legend(g,loc='upper left',
    bbox_to_anchor=[1,1],
    borderaxespad=0,
    handletextpad=.5,
    markerscale=1,
)
#plt.savefig('../plots/pca.pdf',bbox_inches='tight',transparent=True,dpi=300)
# %%
b = meta['mbal']
meta['M']=0
meta.loc[b>=1,'M']=1
meta.loc[b>4,'M']=2
meta.groupby('M').describe()['Death']

# %%
