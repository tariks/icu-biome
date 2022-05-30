# %%
#%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
from scipy import stats

# %%
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)

# %%


pplt.rc.fontsize = 7
pplt.rc["tick.labelsize"] = 7
pplt.rc["grid.linestyle"] = ":"
pplt.rc["label.size"] = 7
pplt.rc["font.small"] = 7
pplt.rc["font.large"] = 7
pplt.rc["abc.size"] = "7pt"
pplt.rc["meta.color"] = "#121212"
pplt.rc["cmap.robust"] = True
pplt.rc["title.border"] = False
pplt.rc["abc.border"] = False
# pplt.rc['cmap.robust'] = False
pplt.rc.textcolor = "#121212"
pplt.rc["meta.color"] = "#121212"
pplt.rc.cycle = "Set1"
pplt.rc.titlepad = 3
pplt.rc["meta.width"] = 0.6
pplt.rc["axes.facecolor"] = "#ffffff"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True
pplt.rc["formatter.zerotrim"] = False
pplt.rc["subplots.align"] = True
pplt.rc["cmap.sequential"] = "Batlow"
pplt.rc["cmap.sequential"] = "Balance"
pplt.rc["colorbar.width"] = 0.08
# pplt.rc.fontfamily = "TeX Gyre Heros"
pplt.rc.fontfamily = "Source Sans Pro"

pplt.rc.inlineformat = "retina"
pplt.config_inline_backend()
blue = (0.000, 0.267, 0.533, 1.000)
red = (0.733, 0.333, 0.400, 1.000)
# %%
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

# %%
def jitter(n, bins=16):
    rang = n.max() - n.min()
    x = pd.cut(
        n,
        bins,
    )
    y = x.copy()
    y = pd.Series(np.zeros(n.size), index=x.index)
    counted = []
    for i in range(n.size):
        count = (x == x[i]).sum()
        if x[i] in counted:
            pass
        elif count == 1:
            pass
        else:
            counted.append(x[i])
            theseguys = y[x == x[i]].index
            width = count * 0.1
            y[theseguys] += np.arange(-width / 2, width / 2, width / theseguys.size)
    y += (np.random.random(y.size) - 0.5) * 0.015
    return y


jitter(meta["Age"])

# %%
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
for k, v in pvals.items():
    print("{}\t{:.3f}".format(k, v))

# %%
fig, axs = pplt.subplots(
    suptitle="Alpha diversity metrics",
    nrows=2,
    ncols=4,
    share=4,
    wspace=0,
    journal="pnas2",
    yloc="left",
    ylabel="Geometric Z-score",
    xloc="bottom",
    gridalpha=0.75,
    gridlinewidth=0.35,
    xminorlocator="none",
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
    s=9,
    mec="#ffffff",
    alpha=0.5,
    marker="x",
    lw=0.5,
)
for metric, ax in zip(div_metrics, axs):
    tmp = meta[metric].copy()
    tmp[:] = stats.gzscore(tmp)
    box1 = (tmp.loc[meta["M"] == 0],)
    box2 = (tmp.loc[meta["M"] == 1],)
    ax.box(
        0,
        tmp.loc[meta["M"] == 0],
        ec=blue,
        **param,
    )
    xcoords = (np.random.random(tmp.loc[meta["M"] == 0].size) - 0.5) * 0.5
    xcoords = jitter(tmp.loc[meta["M"] == 0])
    ax.scatter(x=xcoords, y=tmp.loc[meta["M"] == 0], c=blue, **scatterparam)
    ax.box(
        1,
        tmp.loc[meta["M"] == 1],
        ec=red,
        **param,
    )
    xcoords = 1 + (np.random.random(tmp.loc[meta["M"] == 1].size) - 0.5) * 0.3
    xcoords = 1 + jitter(tmp.loc[meta["M"] == 1])
    ax.scatter(x=xcoords, y=tmp.loc[meta["M"] == 1], c=red, **scatterparam)
    ax.text(
        -0.375,
        -2.2,
        "P={:.3g}".format(pvals[metric]),
        va="top",
    )
    ax.format(
        ylim=(-3, 3), upperlefttitle="  " + metric, xlabel="", xticks="none", titlepad=0
    )
for ax in axs[:, 1:]:
    ax.format(
        ticklen=0,
    )
ax = axs[-1, -1]
ax.format(
    xloc="none",
    yloc="none",
    grid=False,
)
ax.text(0.25, 2, "MMI < 0", va="center")
ax.text(0.25, 1, "MMI > 0", va="center")
ax.hlines([2, 1], [-0.3, -0.3], [0.2, 0.2], c=[blue, red], lw=2)
fig.set(figheight=2.65, figwidth=3.9)

# %%
fig.savefig(
    "../plots/alphabox_ttest.pdf", dpi=600, transparent=True, bbox_inches="tight"
)
