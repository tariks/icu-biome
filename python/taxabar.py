# %%
%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
# %%
pplt.rc.load('../plots/.proplotrc')
pplt.config_inline_backend('retina')
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
mypal = [
'#b6dbff', # lighter blue
'#6db6ff', # light blue
'#009292', # teal
'#006edb', # blue
#'#004949', # dark teal
'#db6d00', # orange
'#924900', # brown
'#920000', # blood red
'#ffff6d', # soft yellow
'#ffb6db', # salmon
'#ff6db6', # hot pink
'#b66dff', # light purple
'#490092', # violet
]
#mypal = sns.color_palette(mypal,desat=1)
#mypal
# %%
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
    lw=0.3,
    #        ec='#292929',
    snap=True,
    aa=True,
    #        alpha=.9,
)
# %%
h = np.zeros(52)
mthres = 0.0

low = [i for i in genus.index if meta.loc[i, "mmi"] < mthres]
high = [i for i in genus.index if meta.loc[i, "mmi"] > mthres]
genus = genus.loc[low + high]
x = np.arange(1, len(low) + 1).tolist() + np.arange(len(low) + 3, 55).tolist()
# genus=genus.loc[meta['mbal2'].sort_values().index]

# %%
z=-.02


fig, ax = pplt.subplot(
    xlim=(0.5, 51.5),
    #figwidth="17.4cm",
    figwidth="8.4cm",
    right=None,
    left=None,
    xticks="none",
    xminorticks="none",
    yminorticks="none",
    xtickloc="none",
    xloc=('data',z),
    refaspect=1.75,
    yloc=('data',-.25)
)

ax.bar(
    x=np.arange(52)+1,
    height=genus.values,
    width=1,
    stack=True,
    cycle=pplt.Cycle(mypal[::-1]),
    ec="#c9c9c9",
    lw=.3,
    labels=genus.columns,
    edgefix=True,
    snap=True,
    alpha=.9,
)
ax.format(
    fontsize=5,
    xloc="none",
    xticks=[],
    yloc="left",
    yformatter="{x:.0%}",
    ylim=(0,1),
)
hparam=dict(
    lw=2,
    solid_capstyle="butt",
    clip_on=False,
    snap=True,
)
ax.axhline(
    y=z,
    xmin=0,
    xmax=21.5/52,
    color="#1770B8",
    **hparam
)
ax.axhline(
    y=z,
    xmin=21.5/52,
    xmax=1,
    color="#B81720",
    **hparam,
)

ax.text(x=11, y=z-0.03, s="MMI < 0", ha="center", va="top", transform="data")
ax.text(x=52.5-15.5, y=z-0.03, s="MMI > 0", ha="center", va="top", transform="data")
ax.margins(0)
h, l = ax.get_legend_handles_labels()
h, l = h[::-1], l[::-1]
leg = ax.legend(h, l, loc="lower left", bbox_to_anchor=(52,0), bbox_transform=ax.transData, ncols=1,
        fontsize=5,frameon=False)
for i in leg.get_patches():
    i.set(lw=0)

print(fig.get_figheight(), fig.get_figwidth())

# %%



# plt.savefig('../plots/taxabars80.pdf',dpi=300,transparent=True,bbox_inches='tight')
fig.savefig("../plots/taxabars_87mm.pdf", dpi=1200, transparent=True,)


# %%
genus.shape
# %%
genus.index=np.arange(52)+1
genus.to_csv('../plots/taxabardata.csv')
# %%
