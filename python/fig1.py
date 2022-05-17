# %%
%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
from skbio.stats.composition import clr, multiplicative_replacement
from scipy.cluster.hierarchy import leaves_list, linkage
from sklearn.decomposition import PCA
from adjustText import adjust_text
from scipy import stats

pplt.rc.fontsize = 7
pplt.rc["tick.labelsize"] = 7
pplt.rc["grid.linestyle"] = ":"
pplt.rc["label.size"] = 7
pplt.rc["font.small"] = 7
pplt.rc["font.large"] = 7
pplt.rc["abc.size"] = '10pt'
pplt.rc["meta.color"] = '#121212'
pplt.rc["cmap.robust"] = True
pplt.rc["title.border"] = False
pplt.rc["abc.border"] = False
# pplt.rc['cmap.robust'] = False
pplt.rc.textcolor = "#121212"
pplt.rc['meta.color'] = "#121212"
pplt.rc.cycle = "Set1"
pplt.rc.titlepad = 3
pplt.rc.inlineformat = 'retina'
pplt.rc["meta.width"] = 0.6
pplt.rc["axes.facecolor"] = "#ffffff"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True
pplt.rc["formatter.zerotrim"] = False
pplt.rc["subplots.align"] = True
pplt.rc["cmap.sequential"] = 'Batlow'
pplt.rc["cmap.sequential"] = 'Balance'
pplt.rc["colorbar.width"] = .08
pplt.config_inline_backend()
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
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)

div_metrics = [
# 'ACE (R)',
 'Chao1 (R)',
 'Fisher (D)',
 'Shannon (D)',
 'Pielou (E)',
 'Dominance',
 'Simpson (D)',
 'Simpson (E)',
]
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
arr = np.array([[1,1,3],[2,2,4],[2,2,0],[5,6,7],[8,9,10]])
arr=np.array([
    [1,3,3],
    [2,3,3],
    [4,4,0],
    [5,6,7],
])
fig,axs = pplt.subplots(arr,share=False,abcloc='l',
                        #journal='pnas2',
                        figwidth=5.5,
                        abctitlepad=5,
                        #refaspect=1.6,
                        #refwidth=None,
                        #refheight=
                        ref=1,
                        hspace=(0,None,None),
                        abc=True,
                        wratios=(1,1,1),
                        hratios=(1.0,1,1.2,.5,),
                        alignx=False,
                        aligny=False,
                        )

## pca plot

x = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
s = (x > 0).sum() / 52
s = s[s.sort_values(ascending=False).index]
s = s[s > 0.2]
x = pd.read_csv("../feature_tables/genus_t80_nozero.csv", index_col=0)
x=x[s.index.to_list()[:30]]
x[:]=multiplicative_replacement(x)
x[:]=clr(x)

pca=PCA(whiten=True,n_components=2)
pca.fit(x)
xproj = pca.transform(x)
loadings=pd.DataFrame(pca.components_.T,index=x.columns,columns=['PC1','PC2'])
l1,l2=loadings['PC1'].copy(),loadings['PC2'].copy()
l1=l1[l1.abs().sort_values(ascending=False).index]
l2=l2[l2.abs().sort_values(ascending=False).index]
L=loadings**2
L=L.sum(axis=1).sort_values(ascending=False)

ax=axs[0]
ax.margins(.05)
ax.scatter(x=xproj[:,0],y=xproj[:,1],
            c=(meta['mmi']>0).astype(int),
            cmap='vlag',
            s=20,
            robust= False,
            #colorbar='r',
            #colorbar_kw={'width': .08, 'length': .9,
            #        'label': 'MMI',
            #        'labelloc': 'right',
            #        'pad': .5,
            #        'lw': 0,
            #        },
            alpha=.7,
)
ax.format(title='PCA biplot',xtickdir='in',xticklabels='none',
)
ax1=axs[1]
ax1.format(ylim=ax.get_ylim(),xlim=ax.get_xlim(),
    xlabel='PC1 (explained variance ={:.2%})'.format(pca.explained_variance_ratio_[0]),
)
for a in [ax,ax1]:
    a.format(
        gridalpha=.7,
        gridlinewidth=.3,
        abcloc='ul',
        #ylabel='PC2 ({:.2%})'.format(pca.explained_variance_ratio_[1]),
        yloc='left',
        xloc='bottom',
    )
taxa=L[:7].index.to_list()
taxa+=[i for i in l1[:2].index.to_list() if i not in taxa]
taxa+=[i for i in l2[:2].index.to_list() if i not in taxa]
texts=[]
for i in taxa:
    a = l1[i]*4
    b = l2[i]*4
    arrows=[]
    arrows.append(ax1.annotate("", xy=(a, b), xytext=(0, 0),
             arrowprops=dict(
                arrowstyle="->",
                shrinkA=0,
                shrinkB=6,
                ),
             ))
    rho,phi=cart2pol(a,b)
    a,b,ha,va=pol2cart(rho+.2,phi)
    texts.append(ax1.text(x=a,y=b,s=i,fontsize=7,ha=ha,va=va))
adjust_text(texts,ax=ax1,autoalign='xy') #,autoalign='xy')
xproj.shape




## -------------------------------------- ##

heatmapdata = pd.read_csv('../plots/heatmap.csv',index_col=0)

ax=axs[2]
ax.format(
    tickminor=False,
    grid=False,
    xloc= ('axes',-.01),
    yloc="right",
    xrotation=90,
    linewidth=0,
    ltitle='Numerator',
    rtitle='Denominator',
    yticklabelsize=6,
)

cols = [
    "Enterobacteriaceae",
    "Anaerococcus",
    "MMI",
    "Parasutterella",
    "Campylobacter",
]

cbk = {
    "lw": 0,
    "ticklabelsize": 7,
    "ticklabels": "{x: .1f}",
    "space": 0.5,
    "tickminor": False,
    'length': .75,
}

ax.heatmap(
    heatmapdata,
    colorbar="l",
    colorbar_kw=cbk,
    cmap="Balance",
    aspect="auto",
    snap=True,
)
ax.hlines(y=np.arange(-.5,31,1),x1=-.5,x2=4.5,lw=.5,ec='#ffffff',snap=True,alpha=.75)
ax.vlines(x=[1.5,2.5],y1=-.5,y2=30.5,lw=1,color='w',snap=True)


## ----------------------------- ##


# taxabars


genus = pd.read_csv('../plots/taxabardata.csv',index_col=0)


ax=axs[3]

ax.bar(
    x=np.arange(52)+1,
    height=genus.values,
    width=1,
    stack=True,
    # cycle='hawaii',
    cycle="twilight_r",
    cycle_kw={"N": 13, 'shift': -1},
    ec="#101010",
    lw=.25,
    labels=genus.columns,
    edgefix=True,
    snap=True,
    alpha=.9,
)
ax.format(
    xloc="none",
    yloc="left",
    yformatter="{x:.0%}",
    linewidth=0,
    ytickwidth=0.5,
    yticklen=0.3,    
    xlim=(0.5, 51.5),
    tickminor=False,
    xticks="none",
    xtickloc="none",
    title='Bacterial compositions on day 1 in the ICU',
)

ax.axhline(
    y=0,
    xmin=0,
    xmax=21.5/52,
    lw=2,
    color="#176AD4",
    solid_capstyle="butt",
    clip_on=False,
)
ax.axhline(
    y=0,
    xmin=21.5/52,
    xmax=1,
    lw=2,
    color="#B71D14",
    solid_capstyle="butt",
    clip_on=False,
)

ax.text(x=11, y=-0.03, s="MMI < 0", ha="center", va="top", transform="data")
ax.text(x=52.5-15.5, y=-0.03, s="MMI > 0", ha="center", va="top", transform="data")

h, l = ax.get_legend_handles_labels()
h, l = h[::-1], l[::-1]
leg = ax.legend(h, l, loc="right", ncols=1, frame=False,
    borderpad=0,
    fontsize=6,
    space=.75,
    pad=0,
    )
print(fig.get_figheight(), fig.get_figwidth())



fig.format(suptitle='Fig1',)
## ----------------------------- ##




param = dict(
    fill=False,
    lw=.5,
)
idx = (meta["mmi"]>0)

for ax,div in zip(axs[4:7],div_metrics[:3]):
    d = {'y': [meta.loc[~idx,div], meta.loc[idx,div]],
         'x': np.array([1,2])}
    ax.boxplot(
        y=meta[div],
        **param
    )
    ax.format(
        xloc='bottom',
        xticks='none',
        xticklabels='none',
        yloc='left',
        ytickdir='in',
        yticklabeldir='in',
        abc=False,
        title=div,
        xlabel='P={:.2g}'.format(stats.ttest_ind(meta.loc[idx,div],meta.loc[~idx,div])[1]),
    )
# %%
