# %%
%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
from skbio.stats.composition import clr,multiplicative_replacement
from scipy.cluster.hierarchy import leaves_list, linkage

# %%
pplt.rc.load('../plots/.proplotrc')
pplt.config_inline_backend('retina')
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
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
red='#490092' # violet
blue='#924900' # brown

# %%


x = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
s = (x > 0).sum() / 52
s = s[s.sort_values(ascending=False).index]
s = s[s > 0.2]
print(s.size)
x = pd.read_csv("../feature_tables/genus_t80_nozero.csv", index_col=0)
balance = ["Enterobacteriaceae", "Anaerococcus", "Parasutterella", "Campylobacter"]
feats = list(set(balance + s.index.to_list()[:30]))
x[:] = multiplicative_replacement(x)
x = x[feats]
x[:] = clr(x)
#x[:] = np.log2(x)
cols=[]
for i in x.columns:
    s = i.replace('[','').replace(']','')
    cols.append(s)
x.columns=cols




# %%
fig, ax = pplt.subplot(
    tickminor=False,
    grid=False,
    yloc= ('axes',-.01),
    xloc="none",
    linewidth=0,
    xticklabels=[],
    left=None,
    right=None,
    figwidth='8.4cm',
    refaspect=4.75,
    ylabelpad=.15,
#    figheight=2.6,
)
cols = [
    "Anaerococcus",
    "Enterobacteriaceae",
    "MMI",
    "Parasutterella",
    "Campylobacter",
]
x["MMI"] = meta["mmi"].copy()
corr = x.corr(method="spearman").loc[x.columns[:-1], cols]

cbk = {
    "ticklabels": "{x: .1f}",
    "space": '.2em',
    "length": .7,
    #'width': .1,
    'c': 'w',
    'align': 'right',
    'lw': .75,
}
for i in cols[:2]+cols[-2:]:
    corr.loc[i,i]=0
corrmax=corr.max().max()
for i in cols[:2]+cols[-2:]:
    corr.loc[i,i]=corrmax

link = linkage(
    corr,
    optimal_ordering=True,
    method='ward',
)
norm = pplt.Norm('diverging', fair=False, vcenter=0)
order = leaves_list(link)

heat = ax.heatmap(
    corr.iloc[order][cols[::-1]].T,
    colorbar="b",
    colorbar_kw=cbk,
    cmap="vlag",
    aspect="auto",
    snap=True,
    norm=norm,
    levels=9,
#    symmetric=True,
)
ax.vlines(x=np.arange(-.5,31,1),y1=-.5,y2=4.5,lw=.5,ec='#ffffff',snap=True,alpha=1)
ax.hlines(y=[-.5,1.5,2.5,4.5],x1=-.5,x2=30.5,lw=3,color='#ffffff',snap=True)
ax.hlines(y=[1.5,2.5],x1=-.47,x2=30.47,lw=.5,color='#000000',snap=True,capstyle='butt',ls=(0,(5,5)))
ax.hlines(y=[.5,3.5],x1=-.5,x2=30.5,lw=.75,color='#ffffff',snap=True,alpha=1)
labs = [i.replace('[','').replace(']','').replace(' group','') for i in corr.index[order]]
for j,i in enumerate(labs):
    c='#121212'
    if i in aerobes:
        c='blueviolet'
    if i in gpac:
        c='orangered'
    ax.text(j-.35,4.42,i,c=c,va='bottom',ha='left',rotation=75,fontsize=5,fontstyle='italic')
#ax.text(.5,30.75,'Numerator',ha='center',va='bottom',fontsize=6.5)
#ax.text(3.5,30.75,'Denominator',ha='center',va='bottom',fontsize=6.5)

#fig.set(figheight=2.6,figwidth='14cm')

#ax.invert_yaxis()
ax.format(xrotation=0,ylabelpad=0)
print(fig.get_figheight(),fig.get_figwidth())
fig.savefig(
    "../plots/taxaheatmap_87mm.pdf",
    dpi=1200,
    transparent=True,
)


# %%
