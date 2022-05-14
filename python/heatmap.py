# %%
#%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
from skbio.stats.composition import clr, multiplicative_replacement
from scipy.cluster.hierarchy import leaves_list, linkage

pplt.rc.fontsize = 7
pplt.rc["tick.labelsize"] = 7
pplt.rc["grid.linestyle"] = ":"
pplt.rc["label.size"] = 7
pplt.rc["cmap.robust"] = True
# pplt.rc['cmap.robust'] = False
pplt.rc.textcolor = "#121212"
pplt.rc.cycle = "Set1"
pplt.rc.titlepad = 9
pplt.rc["meta.width"] = 0.6
pplt.rc["axes.facecolor"] = "#ffffff"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True
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
pplt.rc.fontfamily = "TeX Gyre Heros"
# pplt.rc.fontfamily = "Source Sans Pro"
# %%
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
table = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)

# %%


s = (table > 0).sum() / 52
s = s[s.sort_values(ascending=False).index]
s = s[s > 0.2]
print(s.size)
balance = ["Enterobacteriaceae", "Anaerococcus", "Parasutterella", "Campylobacter"]
feats = list(set(balance + s.index.to_list()[:30]))
# table = pd.read_csv("../feature_tables/genus_t80_nozero.csv", index_col=0)
x = table.copy().loc[:, feats]
# np.random.RandomState(seed=1)
# x[x==0]+=np.random.random(x[x==0].shape)*.9 +.1
x[:] = multiplicative_replacement(x)
# x[:]=clr(x)
x[:] = np.log2(x)

cbk = {"width": 0.1, "lw": 0, "ticklabels": "{x: .1f}", "space": 0}
cmk = {"discrete": True, "N": 7}


# %%
fig, ax = pplt.subplot(
    # journal="aaas2",
    journal="pnas2",
    fontsize=7,
    labelsize=7,
    yticklabelsize=7,
    tickminor=False,
    grid=False,
    xloc="bottom",
    yloc="right",
    xrotation=90,
    linewidth=0,
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
    "width": 0.1,
    "lw": 0.25,
    "ticklabelsize": 7,
    "ticklabels": "{x: .1f}",
    "space": 0.75,
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
    # cmap_kw=cmk,
)
# ax.format(xloc='none',)
fig.savefig(
    "../plots/taxaheatmap_pnas2_helvetica.pdf",
    dpi=600,
    transparent=True,
)

# %%
