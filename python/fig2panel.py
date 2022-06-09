# %%
#%matplotlib inline
import numpy as np
import pandas as pd
from scipy import stats

# %%
import proplot as pplt

# %%
pplt.rc.load('../plots/.proplotrc')
pplt.config_inline_backend('retina')
pplt.rc.fontfamily = "Arial"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True

# %%
array = np.array([
    [1,3,2,],
    [4,5,2,],
    [4,5,0,],
    [6,7,0,],
    [6,7,10,],
    [8,9,10,],
    ]
)
array = np.array([
    [1,3,3,4,5,],
    [6,7,7,8,9,],
    [2,2,10,10,11],
    ]
)
fig, axs = pplt.subplots(
    array,
    yloc="left",
    xloc='bottom',
    hratios=[1,1,1.2],
    wratios=[1,.5,.5,1,1],
    sharex=False,
    sharey=3,
    span=False,
    hspace=('.5em',None),
    wspace='.5em',
    gridalpha=0.7,
    gridlinewidth=0.3,
    refaspect=1,
    figwidth='12.5cm',
)


###### diversity plots
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
blue = (0.000, 0.267, 0.533, 1.000)
red = (0.733, 0.333, 0.400, 1.000)
div_metrics = [
    # 'ACE (R)',
    "Chao1 (R)",
    "Fisher (D)",
    "Shannon (D)",
    "Pielou (E)",
    "Dominance",
    "Simpson (D)",
    "Simpson (E)",
]
meta = meta.iloc[1:, :]


def jitter(n, bins=16):
    rang = n.max() - n.min()
    x = pd.cut(
        n,
        bins,
    )
    y = x.copy()
    y = np.zeros(n.size)
    counted = []
    for i in range(n.size):
        count = (x == x[i]).sum()
        if x[i] in counted:
            pass
        elif count == 1:
            pass
        else:
            counted.append(x[i])
            theseguys = (x == x[i])
            width = count * 0.1
            y[theseguys] += np.arange(-width / 2, width / 2, width / theseguys.sum())
    y += (np.random.random(y.size) - 0.5) * 0.015
    return y


jitter(meta["Age"])

v = "M"
a = meta.loc[meta[v] == 0]
b = meta.loc[meta[v] == 1]
pvals = {
    i: stats.ranksums(
        a[i],
        b[i],
    )[1]
    for i in div_metrics
}
pvals = {
    i: stats.ttest_ind(a[i], b[i], alternative="two-sided", equal_var=False)[1]
    for i in div_metrics
}
divax=axs[:2,:]
divax.format(
    xminorlocator="none",
    ylim=(-3,3),
    xticks='none',
    xlabel='',
    ylabel='',
    ticklen=0,
)
param = dict(
    fill=False,
    means=False,
    meanc="black",
    lw=1.1,
    showcaps=False,
    showfliers=False,
    widths=0.75,
)
scatterparam = dict(
    s=15,
    mec="#ffffff",
    alpha=0.7,
    marker="x",
    lw=0.5,
)
b=[]
for metric, ax in zip(div_metrics, divax[:-1]):
    tmp = meta[metric].copy()
    tmp[:] = stats.gzscore(tmp)
    box1 = (tmp.loc[meta["M"] == 0],)
    box2 = (tmp.loc[meta["M"] == 1],)
    b.append(ax.box(
        0,
        tmp.loc[meta["M"] == 0].values,
        ec=blue,
        autoformat=False,
        label='MMI < 0',
        **param,
    ))
    xcoords = (np.random.random(tmp.loc[meta["M"] == 0].size) - 0.5) * 0.5
    xcoords = jitter(tmp.loc[meta["M"] == 0])
    ax.scatter(x=xcoords, y=tmp.loc[meta["M"] == 0].values, autoformat=False, c=blue, **scatterparam)
    b.append(ax.box(
        1,
        tmp.loc[meta["M"] == 1].values,
        ec=red,
        autoformat=False,
        label='MMI > 0',
        **param,
    ))
    xcoords = 1 + (np.random.random(tmp.loc[meta["M"] == 1].size) - 0.5) * 0.3
    xcoords = 1 + jitter(tmp.loc[meta["M"] == 1].values)
    ax.scatter(x=xcoords, y=tmp.loc[meta["M"] == 1].values, autoformat=False, c=red, **scatterparam)
    ax.text(
        -0.375,
        -2.95,
        "P={:.3g}".format(pvals[metric]),
        va="bottom",
    )
    ax.text(
        -0.375,
        2.15,
        metric,
        va="bottom",
    )
    ax.format(xticks=[])
#subax.format(
#        yticklabels='none',xlabel='',
#    )
divax[0].format(abc=True)
ax = divax[-1]
ax.format(
    xloc="none",
    yloc="left",
    grid=False,
)
bb=[i.get('boxes')[0] for i in b[:2]]
ax.legend(bb,['MMI < 0', 'MMI > 0'],loc='fill',frameon=False,ncols=1,borderpad=1)
ax.spines['left'].set(color=(1,1,1,0))

#ax.text(.35, .6, "MMI < 0", ha="center",transform=ax.transAxes)
#ax.text(.35, .4, "MMI > 0", ha="center",transform=ax.transAxes)
#ax.hlines([.6, .4], [0.05, 0.05], [.3, .3], c=[blue, red], lw=2.5,transform=ax.transAxes)

#subax.format(xticklen=0,xlabel='')
####### pca

from skbio.stats.composition import clr
blue = (0.000, 0.267, 0.533, 1.000)
red = (0.733, 0.333, 0.400, 1.000)
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi
def pol2cart(rho, phi):
    a = np.cos(phi) * rho
    b = np.sin(phi) * rho
    ha = "right" if a < 0 else "left"
    va = "top" if b < 0 else "bottom"
    # if abs(a) > abs(b):
    #    va = "center"
    # else:
    #    ha = "center"
    return a, b, ha, va


meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
mbal = meta["mmi"]

from sklearn.decomposition import SparsePCA

x = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
x.columns = [i.replace('[','').replace(']','') for i in x.columns]
s = (x > 0).sum() / 52
s = s[s.sort_values(ascending=False).index]
s = s[s > 0.2]
x = pd.read_csv("../feature_tables/genus_t80_nozero.csv", index_col=0)
x.columns = [i.replace('[','').replace(']','') for i in x.columns]
# x/=x.sum(axis=0)
# x[:]=multiplicative_replacement(x)
x = x[s.index.to_list()[:30]]
x[:] = clr(x)
# x-=x.mean()

pca = SparsePCA(
    n_components=2,
    random_state=0,
    alpha=1,
    ridge_alpha=0.02,
    max_iter=200,
    method="cd",
)
W = pca.fit_transform(x)
xproj = W
loadings = pd.DataFrame(
    pca.components_.T[:, :2], index=x.columns, columns=["PC1", "PC2"]
)
l1, l2 = loadings["PC1"].copy(), loadings["PC2"].copy()
l1 = l1[l1.abs().sort_values(ascending=False).index]
l2 = l2[l2.abs().sort_values(ascending=False).index]
L = loadings**2
L = L.sum(axis=1).sort_values(ascending=False)
a = meta.loc[meta["month"] == 0].index
b = meta.loc[meta["month"] == 1].index
loadings
subax=axs[-1,:]


handles=[]

subax.format(
    yminorlocator="none",
    ylabel='MMI',
)
W = pd.DataFrame(W, index=x.index, columns=loadings.columns)
l1 = l1.sort_values(ascending=False)
l2 = l2.sort_values(ascending=False)
z = np.zeros(W.shape[0])
cycle=pplt.Cycle('colorblind',ms=10,mew=.5,  )
colormap='batlow_r'
cmap=pplt.Colormap(colormap,left=.3,right=.7)
for ax,pc in zip(subax[:2],['PC1','PC2']):
#    norm=pplt.Normalize(W[pc])
    for i in [0,1]:
        xw=W[pc].loc[meta['month']==i]
        yw=mbal.loc[meta['month']==i]
        if i==1:
            lab='Death < 28 days'
            m='x'
            mew=1.2
            a=1
        else:
            lab='Death > 28 days'
            m='o'
            mew=0
            a=.6

        print(xw.shape)
        handles.append(ax.scatter(
            xw,
            yw,
            cycle=cycle,
            c = (xw-W[pc].min())/W[pc].max(),
            cmap=cmap,
    #        fc=(0,0,0,0),
            marker=m,
            alpha=a,
            mew=mew,
            label=lab,
            snap=True,
        ))
        ax.format(xlabel=pc,) #,ylabelpad=-1.2,xlabelpad=5)
handles=handles[:2]

subax[0].format(abc=True,) #,ylabelpad=-1.2,xlabelpad=5)
#ax.legend(loc='upper center',bbox_to_anchor=[.5,.95+5/13],frameon=False,ncols=1)
blue='dull blue'
red = 'dark peach'
panels=[]
for ax,l in zip(subax[:2],[l1,l2],):
    panel = ax.panel('bottom',share=False,width=.4,yloc='none',ylabel='')
    taxa = l.index.to_list()
    taxa = taxa[:4] + taxa[-4:]
    taxa = taxa[::-1]
    #panel.bar(l.loc[taxa], negpos=True, negcolor=blue, poscolor=red, alpha=0.9, snap=True)
    lo,hi = l.min(),l.max()
    norm=pplt.Normalize()
    c=norm(l.loc[taxa])
    panel.bar(l.loc[taxa], c=cmap(c), snap=True)
    panel.axhline(0, color="#303030", lw=0.8)
    limsx=(-.5,7.5)
    limsy=(-.5,.35)
    panel.format(
        yloc="none",
        xloc="bottom",
        ylabel="",
        xlabel="",
        yminorlocator="none",
        xminorlocator="none",
        yticks=[0],
        grid=False,
        xlim=limsx,
        ylim=limsy,
        xrotation=50,
    )
    panels.append(panel)


print(fig.get_figheight(),fig.get_figwidth())
#axs[4,0].text(-1,0,'Geometric Z-score',rotation=90,ha='right',va='center')
'''
subax.format(
    ylabel='MMI',
)
axs[:,:2].format(
    ylabel='Geometric Z-score',
    xlabel='',
    xticks=[],
)
for ax in axs[:,-1]:
    ax.format(
        ylabel='MMI',
    )
'''   

axs[0].format(ltitle='Diversity metrics in MMI subgroups') 
divax.format(ylabel='Geometric Z-score',yticks=np.arange(-2,3))
axs[-1,0].format(ltitle='MMI vs community structure components (SSPCA)')
panels[0].format(abc=True)
axs[-1].legend(handles,loc='fill',ncols=1,frameon=False,borderpad=1)
fig.savefig('../plots/fig2panel.pdf',dpi=1000,transparent=True)
# %%
