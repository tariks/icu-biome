# %%
%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
from biom import Table
from deicode.rpca import auto_rpca as rpca
from skbio.stats.composition import clr, multiplicative_replacement
from scipy.cluster.hierarchy import leaves_list, linkage

pplt.rc.fontsize = 7
pplt.rc["tick.labelsize"] = 7
pplt.rc["grid.linestyle"] = ":"
pplt.rc["label.size"] = 7
pplt.rc["font.small"] = 7
pplt.rc["font.large"] = 7
pplt.rc["abc.size"] = '7pt'
pplt.rc["meta.color"] = '#121212'
pplt.rc["cmap.robust"] = True
pplt.rc["title.border"] = False
pplt.rc["abc.border"] = False
# pplt.rc['cmap.robust'] = False
pplt.rc.textcolor = "#121212"
pplt.rc['meta.color'] = "#121212"
pplt.rc.cycle = "Set1"
pplt.rc.titlepad = 3
pplt.rc["meta.width"] = 0.6
pplt.rc["axes.facecolor"] = "#ffffff"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True
pplt.rc["formatter.zerotrim"] = False
pplt.rc["subplots.align"] = True
pplt.rc["cmap.sequential"] = 'Batlow'
pplt.rc["cmap.sequential"] = 'Balance'
pplt.rc["colorbar.width"] = .08
gpac = [
    "Anaerococcus",
    "Fenollaria",
    "Finegoldia",
    "Peptococcus",
    "Peptostreptococcus",
    "Coprococcus",
    "Atopobium",
    "Ruminococcus",
    "Parvimonas",
    "Peptoniphilus",
    "Blautia",
    "Gallicola",
    "Murdochiella",
    "Sarcina",
]
aerobes = [
    "Staphylococcus",
    "Enterobacteriaceae",
    "Lactobacillus",
    "Enterococcus",
    "Acinetobacter",
    "Pseudomonas",
    "Streptococcus",
]
#pplt.rc.fontfamily = "TeX Gyre Heros"
pplt.rc.fontfamily = "Source Sans Pro"

pplt.config_inline_backend()
pplt.rc.inlineformat = 'retina'

# %%
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi

def pol2cart(rho, phi):
    a = np.cos(phi) * rho
    b = np.sin(phi) * rho
    ha = "right" if a < 0 else "left"
    va = "top" if b < 0 else "bottom"
    #if abs(a) > abs(b):
    #    va = "center"
    #else:
    #    ha = "center"
    return a, b, ha, va


# %%
table = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
mbal = meta["mmi"]
#norm = plt.cm.colors.CenteredNorm(vcenter=0, halfrange=2, clip=True)
#norm = plt.cm.colors.Normalize(vmin=-3, vmax=3, clip=True)

# %%
table = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
T = Table(table.T.values, table.columns, table.index)
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
oord, dis = rpca(T, min_feature_frequency=20, min_feature_count=5)
pc, L, ex = oord.samples, oord.features, oord.proportion_explained
meta=meta.loc[meta.index.isin(pc.index)]
top = L["PC1"] ** 2 + L["PC2"] ** 2
top = top.sort_values()[::-1]
print(top[:10], top.shape[0])
a = L["PC1"][L["PC1"].abs().sort_values()[::-1].index]
b = L["PC2"][L["PC2"].abs().sort_values()[::-1].index]
topa, topb = a.index[0], b.index[0]
meta['PC1'] = pc.loc[meta.index,'PC1']
meta['PC2'] = pc.loc[meta.index,'PC2']
meta.shape


# %%
fig, ax = pplt.subplot(
    journal="nat1",
    tickminor=False,
    xloc='bottom',
    yloc='left',
)
ax.scatter(x=meta['PC1'],y=meta['PC2'],
           c=(meta['mmi']>0).astype(int),colorbar='r',
           cmap='batlow',
           alpha=.9,
           snap=True,
           aa=True,
)
ax.format(
    xlabel='PC1 (explained variance ={:.2%})'.format(oord.proportion_explained[0]),
    ylabel='PC2 (explained variance ={:.2%})'.format(oord.proportion_explained[1]),
    gridalpha=.6,
    gridlinewidth=.3,
)

# %%




from sklearn.decomposition import SparsePCA, PCA
from adjustText import adjust_text
x = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
s = (x > 0).sum() / 52
s = s[s.sort_values(ascending=False).index]
s = s[s > 0.2]
x = pd.read_csv("../feature_tables/genus_t80_nozero.csv", index_col=0)
#x/=x.sum(axis=0)
#x[:]=multiplicative_replacement(x)
x=x[s.index.to_list()[:30]]
x[:]=clr(x)

#x[:]=np.log10(x)
#x[x>0]=np.log2(x[x>0])
#x[:]=quantile_transform(x,)
#x[:]=power_transform(x)

#x[:]=StandardScaler().fit_transform(x)
#pca=PCA(whiten=True,n_components=2)
#pca=NMF(
#    n_components=2,
#        init='nndsvdar',
#        init='random',
#        solver='mu',
#        beta_loss = 'kullback-leibler',
#        tol = .0001,
#        max_iter=100,
#        random_state=0,
#        l1_ratio=.5,
#        alpha=.0001,
#        )
pca=SparsePCA(
    n_components=2,
    random_state=0,
    alpha=3,
    ridge_alpha=.01,
    max_iter=200,
    method='cd',

)
W = pca.fit_transform(x)
xproj = W
loadings=pd.DataFrame(pca.components_.T[:,:2],
    index=x.columns,
    columns=['PC1','PC2'])
l1,l2=loadings['PC1'].copy(),loadings['PC2'].copy()
l1=l1[l1.abs().sort_values(ascending=False).index]
l2=l2[l2.abs().sort_values(ascending=False).index]
L=loadings**2
L=L.sum(axis=1).sort_values(ascending=False)
loadings




# %%
fig,axs=pplt.subplots(
    figtitle='Decomp',
    nrows=2,
    ncols=1,
    share=3,
    hspace=0,
    journal='nat1',
    yloc='left',
    xloc='bottom',
    gridalpha=.7,
    gridlinewidth=.3,

    )
ax=axs[0]
ax.margins(.05)
ax.scatter(x=xproj[:,0],y=xproj[:,1],
#            c=(mbal>0).astype(int),
            c=mbal,
            #cmap='vik',
            #cmap='vlag',
            cmap='Batlow',
            robust= False,
            colorbar='right',
            colorbar_kw={'width': .08, 'length': .9,
                    'label': 'MMI',
                    'labelloc': 'right',
                    'pad': .5,
                    },
            alpha=.7,
)
ax.format(xtickdir='in',)
fig.format(
#    xlabel='PC1 (explained variance ={:.2%})'.format(pca.explained_variance_ratio_[0]),
#    ylabel='PC2 (explained variance ={:.2%})'.format(pca.explained_variance_ratio_[1]),
)
ax1=axs[1]
ax1.format(ylim=ax.get_ylim(),xlim=ax.get_xlim(),)
taxa=L[:5].index.to_list()
#taxa+=[i for i in l1[:2].index.to_list() if i not in taxa]
#taxa+=[i for i in l2[:2].index.to_list() if i not in taxa]
texts=[]
arrows=[]
maxcomponent = max(l1.abs().max(),l2.abs().max())
maxcomponent=loadings.mean().mean()
limx = ax.get_xlim()
limx = (limx[0]-.1,limx[1])
limy = ax.get_ylim()
limy = (limy[0]-.1,limy[1])

#fig.format(xlim=limx,ylim=limy)
for i in taxa:
    a = l1[i]
    b = l2[i]
#    a/=limx[1]
#    b/=limy[1]
    arrows.append(ax1.arrow(0,0,a,b,
        width=.005,
        head_width=.01,
        length_includes_head=True,
        color='black',
        alpha=1,
        ))

    rho,phi=cart2pol(a,b)
    a,b,ha,va=pol2cart(rho+.2,phi)
    texts.append(ax1.text(x=a,y=b,s=i,fontsize=7,ha=ha,va=va))
adjust_text(texts,ax=ax1,autoalign='xy') #,autoalign='xy')

# %%
fig.savefig('../plots/pca_genus80.pdf',dpi=600,bbox_inches='tight')
# %%
def roll(sides=20,n=1):
    for i in range(n):
        print(np.random.randint(1,sides+1))

# %%
blue = (0.000, 0.467, 0.733, 1.000)
red = (0.933, 0.467, 0.200, 1.000)

# %%
fig,axs=pplt.subplots(
    suptitle='Principle components vs MMI',
    nrows=2,
    ncols=1,
    sharex=0,
    sharey=False,
    hspace=0,
    journal='nat1',
    yloc='left',
    ylabel='PC1',
    xloc='bottom',
    xlabel='MMI',
    gridalpha=.7,
    gridlinewidth=.3,
    xminorlocator='none',
    yminorlocator='none',

    )
ax=axs[0]
W=pd.DataFrame(W,index=x.index,columns=loadings.columns)
l1=l1.sort_values(ascending=False)
l2=l2.sort_values(ascending=False)
limsx=(loadings.min().min()*1.1,loadings.max().max()*1.1)
limsy=(W.min().min()*1.1,W.max().max()*1.1)
#ax.text(0.05,.98,'\n'.join(l1.index.to_list()[:3]),transform='axes',va='top')
#ax.text(0.05,.02,'\n'.join(l1.index.to_list()[-3:]),transform='axes',va='bottom')
z=np.zeros(W.shape[0])
ax.vlines(mbal,W['PC1'],negpos=True,negcolor=blue,poscolor=red,
    lw=.85,snap=True,
    )

panel = ax.panel('right',
        width=.65,
        share=False,
        )
taxa=l1.index.to_list()
taxa = taxa[:4]+taxa[-4:]
taxa=taxa[::-1]
panel.barh(l1.loc[taxa],negpos=True,negcolor=blue,poscolor=red,alpha=.7,snap=True)
panel.format(yloc='right',xloc='bottom',ylabel='',xlabel='',yminorlocator='none',xminorlocator='none',
            xlim=limsx,xticks=[0])
panel.axvline(0,color='#333333',lw=.8)

ax=axs[1]
#ax.text(0.05,.98,'\n'.join(l2.index.to_list()[:3]),transform='axes',va='top')
#ax.text(0.05,.02,'\n'.join(l2.index.to_list()[-3:]),transform='axes',va='bottom')
#ax.scatter(W[:,1],mbal,c=W[:,0])
z=np.zeros(W.shape[0])
ax.vlines(mbal,W['PC2'],negpos=True,negcolor=blue,poscolor=red,
    lw=.85,snap=True,
    )
ax.format(ylabel='PC2')
panel = ax.panel('right',
        width=.65,
        share=False,
        )
taxa=l2.index.to_list()
taxa = taxa[:4]+taxa[-4:]
taxa=taxa[::-1]
panel.barh(l2.loc[taxa],negpos=True,negcolor=blue,poscolor=red,alpha=.7,snap=True)
panel.format(yloc='right',xloc='bottom',ylabel='',xlabel='coef.',yminorlocator='none',xminorlocator='none',
            xticks=[0],xlim=limsx)
panel.axvline(0,color='#333333',lw=.8)

for i in axs:
    i.format(ylim=limsy)
fig.set(figheight=3,figwidth=3.75)
# %%

fig.savefig('../plots/sparse_pca_genus80.pdf',dpi=600,bbox_inches='tight')
# %%
