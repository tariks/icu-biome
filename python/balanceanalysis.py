# %%
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
from skbio.stats.composition import multiplicative_replacement, clr
#from sklearn.preprocessing import PowerTransformer,QuantileTransformer 
import cmasher as cmr
from pca import pca
from deicode.rpca import auto_rpca as rpca
from biom import Table
plt.style.use("../lance.txt")
gpac=[
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
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi
def pol2cart(rho,phi):
        a=np.cos(phi)*rho
        b=np.sin(phi)*rho
        ha = 'right' if a<0 else 'left'
        va = 'top' if b<0 else 'bottom'
        if abs(a) > abs(b):
            va='center'
        else:
            ha='center'
        return a,b,ha,va
# %%
meta = pd.read_csv("../meta/52_bal.csv", index_col=0)
#table = pd.read_csv("../feature_tables/52_phyla70.csv", index_col=0)
#table = pd.read_csv("../feature_tables/52N_genus_t80_full.csv", index_col=0)
table = pd.read_csv("../feature_tables/52N_genus_t70.csv", index_col=0)
#table = pd.read_csv('../feature_tables/all_asv.csv',index_col=0)
#table=table.loc[table.index.isin(meta.index)]
T = Table(table.T.values,table.columns,table.index)
mbal=meta['mbal']
#norm=plt.cm.colors.CenteredNorm(vcenter=2.1,halfrange=4)
norm=plt.cm.colors.Normalize(vmin=-4,vmax=6.5,clip=True)
#cmap=plt.cm.get_cmap('coolwarm')
cmap=cmr.dusk
#cmap=cmr.ember
arrow=dict(
#    head_width=.3,
#    tail_width=.6,
    width=.5,
    headwidth=3,
    headlength=3,
    color='#810004',
    alpha=.5,
)
# %%
oord,dis = rpca(T,min_feature_frequency=10,min_feature_count=200)
pc,L,ex = oord.samples, oord.features, oord.proportion_explained
top=L['PC1']**2 + L['PC2']**2
top= top.sort_values()[::-1]
print(top[:10],top.shape[0])
a = L['PC1'][L['PC1'].abs().sort_values()[::-1].index]
b = L['PC2'][L['PC2'].abs().sort_values()[::-1].index]
c=cmap(norm(mbal[pc.index]))
#c=plt.cm.get_cmap('coolwarm')(norm(mbal[pc.index]))
topa,topb=a.index[0],b.index[0]
fig,ax=plt.subplots()
ax.scatter(pc['PC1'],pc['PC2'],
    c=c,
    marker='o',s=20,
    lw=0,snap=True,aa=True,alpha=.8,
)
#for tax in [topa,topb]:
for tax in top.index[:3]:
    x,y=.7*a[tax],.7*b[tax]
    ax.annotate('',(x,y),xytext=(0,0),
    arrowprops=arrow,)
    rho,phi=cart2pol(x,y)
    x,y,ha,va=pol2cart(rho+.02,phi)
    ax.text(x,y,tax,
        ha=ha,va=va,fontsize=5)

ax.set_xlabel('PC1, expl. var. ={:.2%}'.format(ex[0]),loc='right')
ax.set_ylabel('PC2, expl. var. ={:.2%}'.format(ex[1]),loc='top')
ax.set_title('Robust Aitchison PCA',loc='left')
ax.margins(.05)
limx,limy = ax.get_xlim(),ax.get_ylim()
ax.set_ylim(min(limy[0],.6*b[topb]-.1),max(limy[1],.6*b[topb]+.1))
ax.set_xlim(min(limx[0],.6*a[topa]-.2),max(limx[1],.6*a[topa]+.2))
cbar=plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax,
    fraction=.03,location='right',aspect=25,shrink=.6,pad=.01
    )
cbar.outline.set(lw=0)
plt.savefig('../plots/pca_genus.pdf',dpi=300,transparent=True,bbox_inches='tight')
# %%
fig,axs=plt.subplots(2,1,sharex=False)
fig.set(figheight=6,figwidth=2)
fig.subplots_adjust(hspace=.2)
a=a[:20].sort_values()[::-1]
b=b[:20].sort_values()[::-1]
ax.margins(.05)
ax=axs[0]
ax.set_xlabel('Magnitude of projection on PC1')
ax.set_title('Top 20 taxa by PC loading')
ax.barh(y=a.index,width=a)
ax=axs[1]
ax.barh(y=b.index,width=b)
ax.set_xlabel('Magnitude of projection on PC2')
ax.margins(.05)
plt.savefig('../plots/pca_loadings.pdf',dpi=300,transparent=True,bbox_inches='tight')
# %%
#x=table.copy()
#x[:]=multiplicative_replacement(x)
#x[:]=clr(x)
#x=pd.DataFrame(ilr(x),index=x.index)
#top=x.mean().sort_values()[::-1].index[:].to_list()
#x[:]=PowerTransformer(standardize=True).fit_transform(x)
#x[:]=QuantileTransformer().fit_transform(x)
#x[:]=StandardScaler().fit_transform(x)
x=rpca(T)
# %%
mod = pca(normalize=True)
results = mod.fit_transform(x[top])
fig,ax = mod.biplot(
#    y = (meta['mbal']>1.9).astype(int),
    n_feat=4,
    PC=[1,2],
#    SPE=True,
#    hotellingt2=True,
    figsize=(2.5,2.5),
#    alpha_transparency=.9,
    visible=False,
    color_arrow='#810004',
    verbose=0,
    title='PCA of top 100 taxa',
    legend=False,
    label=None,
)
for t in ax.texts:
    if len(t.get_text())<5:
        t.remove()
    if '(' in t.get_text():
        s = t.get_text().split('(')
        number = float(s[1].split(')')[0])
        s=s[0]+'({:.2f})'.format(number)
        a,b=t.get_position()
        rho,phi=cart2pol(a/1.11,b/1.11)
        rho+=.2
        a=np.cos(phi)*rho
        b=np.sin(phi)*rho
        ha = 'right' if a<0 else 'left'
        va = 'top' if b<0 else 'bottom'
        if abs(a) > abs(b):
            va='center'
        else:
            ha='center'
        plt.setp(t,text=s,color='#333333',fontsize=5,x=a,y=b,ha=ha,va=va)
for i,marker in enumerate(ax.collections):
    marker.set_sizes(marker.get_sizes()/6)
#    marker.set(color='w',lw=.8,
    marker.set(lw=0, color=colors,)
#ax.collections[1].set_edgecolor('#810004')
#ax.collections[0].set_edgecolor('#004481')
ax.set_visible(True)
fig.set_visible(True)
ax.plot()
'''
ax.legend(labels=['Low balance','High balance'],
    loc='upper right',
    bbox_to_anchor=[.99,1.125],
    markerscale=.99,
    handletextpad=.1,
)
'''
cbar=plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='coolwarm'), ax=ax,
    fraction=.03,location='right',aspect=25,shrink=.6,pad=.01
    )
cbar.outline.set(lw=0)
#plt.savefig('../plots/pca_phylum.pdf',dpi=300,transparent=True,bbox_inches='tight')
#plt.savefig('../plots/pca_genus.pdf',dpi=300,transparent=True,bbox_inches='tight')
plt.show()
results['topfeat'][:10]

# %%
mod.plot()
# %%
meta['M']=0
meta.loc[meta['mbal']>1.9,'M']=1
#meta.loc[meta['mbal']>4,'M']=2
meta.groupby("M").describe()["Death"]
# %%
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

b = meta['mbal']
meta['M']=0
meta.loc[b>=1,'M']=1
meta.loc[b>4,'M']=2
meta.groupby('M').describe()['Death']

# %%

#y=power_transform(x.iloc[4,:].values.reshape(-1,1),)
y=x.iloc[4,:]
sns.histplot(y)
# %%
l=results['loadings'].loc['PC1']
l[l.abs().sort_values().index[::-1]]

# %%
