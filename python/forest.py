# %%
#%matplotlib inline
import proplot as pplt
import numpy as np
import pandas as pd
import lifelines

# %%
pplt.rc.load('../plots/.proplotrc')
pplt.config_inline_backend('retina')
pplt.rc.fontfamily = "Arial"
pplt.rc["text.antialiased"] = True
pplt.rc["lines.antialiased"] = True
# %%
def sigstar(n):
    x = int(-np.log10(n))
    return "*" * x


meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
x = pd.read_csv("../survival/univariate.csv", index_col=0)
vdict = pd.read_csv("../meta/vdict.csv", index_col=0, header=None)
vdict_r = {vdict.loc[i, 1]: i for i in vdict.index}
vdict = {v: k for k, v in vdict_r.items()}

# %%
kparam = dict(ci_show=False, xlabel="", yticks=[], aa=True)
kparam = dict(
    ci_show=False,
    xlabel="",
    grid=True,
    aa=True,
    lw=1.2,
)


def plotkm2(ax, v="M"):
    a = meta.loc[meta[v] == 1]
    b = meta.loc[meta[v] == 0]
    km = lifelines.KaplanMeierFitter(alpha=0.1)
    km.fit(a["Death"], a["month"])
    km.plot(**kparam, ax=ax, label=vdict[v], )
    km.fit(b["Death"], b["month"])
    km.plot(
        **kparam,
        ax=ax,
        label=vdict[v].replace(">", "<"),
        legend=False,
    )
    for i in [0,10,20,30]:
        asum=(a['Death']>i).sum()
        bsum=(b['Death']>i).sum()
        ax.text(i,.21,str(asum),ha='center',color="#9F1800")
        ax.text(i,.26,str(bsum),ha='center',color="blue")


# %%
x = x.loc[x.index[::-1]]
x.index = [vdict[i] for i in x.index]
idx = ["Ciprofloxacin", "Ceftriaxone", "Piperacillin-tazobactam", "Vancomycin"]
idx += [
    "Diabetes",
    "Sepsis-3",
    "Immunosuppressed",
    "Malignancy",
    "ARDS",
    "APACHE II score",
]
idx += [
    "Gender (male)",
    "Age",
]
idx += ["MMI > 0", "MMI"]
x = x.loc[idx]
x
# %%
yticks = np.arange(0, 20).tolist()
yticks.remove(4)
yticks.remove(5)
yticks.remove(12)
yticks.remove(13)
yticks.remove(16)
yticks.remove(17)

# %%
fig, axs = pplt.subplots(
    nrows=1,
    ncols=2,
    figwidth='13cm',
    share=False,
    fontsize=7,
    labelsize=7,
    ticklabelsize=7,
    gridalpha=.6,
    gridstyle=':'
)
fig.format(
    linewidth=0.6,
    yticks=pplt.Locator("null"),
    yticklen=0,
    xloc="none",
    yloc="left",
)
fig.set(
    figheight=2.75,
)
ax = axs[0, 0]
ax.hlines(
    x1=x["exp(coef) lower 95%"],
    x2=x["exp(coef) upper 95%"],
    y=yticks,
    color="#121212",
    snap=True,
    aa=True,
    lw=1,
    alpha=0.5,
)
ax.axvline(1, color="#121212", snap=True, aa=True, lw=0.6, alpha=0.6)
ax.scatterx(
    x=x["exp(coef)"],
    y=yticks,
    c="#333333",
    s=22,
    marker="D",
    snap=True,
    aa=True,
    alpha=0.5,
    mew=0,
)
ylim = list(ax.get_ylim())
ylim[1] += 1


ax.format(
    xlim=(-0, 4.5),
    grid=False,
    xloc=("axes", -0.06),
    ylim=ylim,
    xlabel="Hazard ratio",
    title="Univariate Cox analysis",
    titleloc="left",
    yloc="none",
)




panel = ax.panel("left", width=1.3, space='2em')

for i, v in enumerate(x.index):
    sigs = ""
    if x.loc[v, "p"] < 0.05:
        sigs = sigstar(x.loc[v, "p"].astype(float))
    tmp = vdict_r[v]
    if meta[tmp].max() == 1:
        panel.text(
            0.25,
            yticks[i],
            s="  (n={}) ".format(meta[tmp].sum().astype(int)) + str(v) + sigs,
            va="center",
            ha="left",
        )  # ,fontsize=7)
    else:
        panel.text(
            0.27, yticks[i], s=str(v) + sigs, va="center", ha="left"
        )  # ,fontsize=7)

panel.format(
    yticklabels=pplt.Formatter("null"),
    yticklen=0,
    xlim=(-0.02, 2),
    ylim=ax.get_ylim(),
    grid=False,
    xloc="none",
    yloc="none",
)
tparam = dict(
    weight="demi",
    ha="left",
    va="center",
)
panel.text(0, 20.3, s="Microbiome variables", **tparam)
panel.text(0, 4.3, s="Antibacterials", **tparam)
panel.text(0, 12.3, s="Clinical variables", **tparam)
panel.text(0, 16.3, s="Demographic variables", **tparam)

ax = axs[0, 1]

ax.format(
    yloc=("axes", -0.08),
    ylim=(0.2, 1),
    xlim=(0, 40),
    yticks="auto",
    yformatter="{x:.0%}",
    fontsize=7,
    labelsize=7,
    ticklabelsize=7,
    xloc=("axes", -0.06),
    xlabel="",
    lw=0.6,
    grid=True,
    title="Kaplan-Meier curves",
    ylabel="Survivorship",
)
ax.twiny(
    # xcolor="crimson",
    xlabel="Days after ICU admission",
    xlim=ax.get_xlim(),
)
plotkm2(ax=ax)
h, l = ax.get_legend_handles_labels()
leg = ax.legend(
    h,
    l,
    loc='upper right',
    ncols=1,
    fontsize=7,
    frameon=False,
)
print(fig.get_figwidth(), fig.get_figheight())
for t in ax.texts:
    if "time" in t.get_text():
        t.remove()
for i in ax.lines:
    i.set(lw=1, snap=True, aa=True)

# drawmedian(ax,)

# %%
fig.savefig("../plots/forest.pdf", dpi=1000, transparent=True, bbox_inches="tight")
# %%
"""
nat1 = 3.5 1.55
nat2 = 7.2 2.79
agu2 = 7.48 4.52
pnas2 = 4.49 1.88
aaas2 = 4.72 1.96
ams3 = 5.5 2.2
ams2 = 4.5 1.89
"""

# %%
def q13(v):
    q1 = int(v.quantile(0.25))
    q3 = int(v.quantile(0.75))
    return "{},{}".format(q1, q3)


def drawmedian(ax, vector=meta["Death"], group="M", where="left"):
    med = vector.median()
    mediqr = q13(vector)
    pos = vector[meta[group] == 1].median()
    posiqr = q13(vector[meta[group] == 1])
    neg = vector[meta[group] == 0].median()
    neqiqr = q13(vector[meta[group] == 0])
    set1 = pplt.get_colors("Set1")[:2]
    colors = [set1[0], "#333333", set1[1]]
    ax.vlines(x=[pos, med, neg], y1=-0.04, y2=0.01, lw=1, alpha=0.8, colors=colors)
    ax.text(
        x=pos,
        y=0,
        s="Median[IQR]={}[{}]".format(pos, posiqr),
        color=colors[0],
        ha="right",
        transform="data",
    )
    ax.text(
        x=med,
        y=0,
        s="Median[IQR]={}[{}]".format(med, mediqr),
        color=colors[1],
        ha="left",
        transform="data",
    )
    ax.text(
        x=neg,
        y=0,
        s="Median[IQR]={}[{}]".format(neg, negiqr),
        color=colors[2],
        ha="left",
        transform="data",
    )


# %%
