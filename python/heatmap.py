# %%
%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
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

# %%


x = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
s = (x > 0).sum() / 52
s = s[s.sort_values(ascending=False).index]
s = s[s > 0.2]
print(s.size)
x = pd.read_csv("../feature_tables/genus_t80_nozero.csv", index_col=0)
balance = ["Enterobacteriaceae", "Anaerococcus", "Parasutterella", "Campylobacter"]
feats = list(set(balance + s.index.to_list()[:30]))
#x[:] = multiplicative_replacement(x)
x = x[feats]
x[:] = clr(x)
#x[:] = np.log2(x)
cols=[]
for i in x.columns:
    s = i.replace('[','')
    s = i.replace(']','')
    cols.append(s)
x.columns=cols

cbk = {"width": 0.1, "lw": 0, "ticklabels": "{x: .1f}", "space": 0}


# %%
fig, ax = pplt.subplot(
    # journal="aaas2",
    journal="pnas2",
    tickminor=False,
    grid=False,
    xloc= ('axes',-.01),
    yloc="right",
    xrotation=90,
    linewidth=0,
    ltitle='Numerator taxa',
    rtitle='Denominator taxa',
    title='MMI',
)
cols = [
    "Enterobacteriaceae",
    "Anaerococcus",
    "MMI",
    "Parasutterella",
    "Campylobacter",
]
x["MMI"] = meta["mmi"].copy()
corr = x.corr(method="spearman").loc[x.columns[:-1], cols]

cbk = {
    "lw": 0,
    "ticklabelsize": 7,
    "ticklabels": "{x: .1f}",
    "space": 0.5,
    "tickminor": False,
}

link = linkage(
    corr,
    optimal_ordering=True,
    # metric='euclidean',
    metric="braycurtis",
    # method='ward',
    method="complete",
)
order = leaves_list(link)
heat = ax.heatmap(
    corr.iloc[order][cols],
    colorbar="l",
    colorbar_kw=cbk,
    cmap="Balance",
    aspect="auto",
    snap=True,
)
ax.hlines(y=np.arange(-.5,31,1),x1=-.5,x2=4.5,lw=.5,ec='#ffffff',snap=True,alpha=.75)
ax.vlines(x=[1.5,2.5],y1=-.5,y2=30.5,lw=1,color='w',snap=True)
fig.savefig(
    "../plots/taxaheatmap_pnas2.pdf",
    dpi=600,
    transparent=True,
    bbox_inches='tight'
)

# %%
corr = corr.iloc[order][cols]
corr

# %%
corr.to_csv('../plots/heatmap.csv')
# %%
