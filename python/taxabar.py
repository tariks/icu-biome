# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import proplot as pplt

pplt.rc.fontsize = 7
pplt.rc["tick.labelsize"] = 7
pplt.rc["grid.linestyle"] = ":"
pplt.rc["label.size"] = 7
# pplt.rc['cmap.robust'] = True
pplt.rc.textcolor = "#121212"
pplt.rc.cycle = "Set1"
pplt.rc.titlepad = 9
pplt.rc["meta.width"] = 0.6
pplt.rc["axes.facecolor"] = "#ffffff"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True

# pplt.rc.fontfamily='TeX Gyre Heros'
pplt.rc.fontfamily = "Source Sans Pro"
# %%
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
mypal = [
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
# mypal = sns.color_palette(mypal,desat=1)
# %%
%matplotlib inline
pplt.config_inline_backend()
# plt.style.use("../lance.txt")
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
taxmap = pd.read_csv("../feature_tables/taxmap.csv", index_col=0, header=None)
phygen = {
    "Bacteroidetes": [
        "Bacteroides",
        "Parabacteroides",
        "Alistipes",
        "Other Bacteroidetes",
    ],
    "Proteobacteria": [
        "Enterobacteriaceae",
        "Pseudomonas",
        "Acinetobacter",
        "Other Proteobacteria",
    ],
    "Firmicutes": ["GPAC", "Enterococcus", "Staphylococcus", "Other Firmicutes"],
}
colors = {
    #    'Firmicutes': sns.color_palette('tab20b',4)[-4:][::-1],
    #    'Proteobacteria': sns.color_palette('tab20b',12)[-4:][::-1],
    #    'Bacteroidetes': sns.color_palette('tab20c',4),
    "Firmicutes": mypal[4:8][::-1],
    "Proteobacteria": mypal[-4:],
    "Bacteroidetes": mypal[:4][::-1],
}
taxa = []
c = []
for i in ["Proteobacteria", "Firmicutes", "Bacteroidetes"]:
    taxa += phygen.get(i)[::-1]
    c += list(colors.get(i))
meta["mmi"].describe()
# %%
genus = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
# genus = pd.read_csv("../feature_tables/52N_genus_t70.csv", index_col=0)
genus["GPAC"] = 0
for i in genus.columns:
    if i in gpac:
        genus["GPAC"] += genus[i]
        del genus[i]
nont = [i for i in genus.columns if i not in taxa]
for i in ["Bacteroidetes", "Proteobacteria", "Firmicutes"]:
    genus["Other " + i] = 0
    for t in nont:
        if taxmap.loc[t].values[0] == i:
            genus["Other " + i] += genus[t]

genus = genus[taxa]
genus = genus.divide(genus.sum(axis=1), axis=0)
allBacter = genus[phygen["Bacteroidetes"]].sum(axis=1)
order = allBacter.sort_values(ascending=False).index
genus = genus.loc[order]

params = dict(
    width=1,
    lw=0.25,
    #        ec='#292929',
    snap=True,
    aa=True,
    #        alpha=.9,
)
# %%
h = np.zeros(52)
mthres = 0.0

# fig.set(figheight=2.5,figwidth=4)
low = [i for i in genus.index if meta.loc[i, "mmi"] < mthres]
high = [i for i in genus.index if meta.loc[i, "mmi"] > mthres]
genus = genus.loc[low + high]
x = np.arange(1, len(low) + 1).tolist() + np.arange(len(low) + 3, 55).tolist()
# genus=genus.loc[meta['mbal2'].sort_values().index]

# %%

fig, ax = pplt.subplot(
    xlim=(0.5, 51.5),
    journal="pnas2",
    fontsize=7,
    labelsize=7,
    ticklabelsize=7,
    tickminor=False,
    xticks="none",
    xtickloc="none",
)


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
    fontsize=7,
    xloc="none",
    yloc="left",
    yformatter="{x:.0%}",
    linewidth=0,
    ytickwidth=0.5,
    yticklen=0.3,
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

# leg=ax.legend(loc='right',ncols=1,fontsize=7,frame=False)
h, l = ax.get_legend_handles_labels()
h, l = h[::-1], l[::-1]
# leg.remove()
leg = ax.legend(h, l, loc="right", ncols=1, fontsize=7, frame=False)
fig.set_figheight(2)
print(fig.get_figheight(), fig.get_figwidth())

# %%



# plt.savefig('../plots/taxabars80.pdf',dpi=300,transparent=True,bbox_inches='tight')
fig.savefig("../plots/taxabars80.pdf", dpi=600, transparent=False, bbox_inches="tight")


# %%
genus.shape
# %%
genus.index=np.arange(52)+1
genus.to_csv('../plots/taxabardata.csv')
# %%
