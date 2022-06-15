# %%
#%matplotlib inline
import proplot as pplt
import numpy as np
import pandas as pd
import lifelines

# %%
pplt.rc.load('../plots/.proplotrc')
pplt.config_inline_backend('retina')
red,blue = '#C91923','#0C385C'
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
kparam = dict(
    lw=.8,
    drawstyle='steps-post',
    solid_capstyle='butt',
    solid_joinstyle='round',
    snap=True,
    autoformat=False,
)


def plotkm2(ax, v, c="#999999",ls='-', label='Cohort (n=52)',offset=.01):
    a = meta.loc[v]
    km = lifelines.KaplanMeierFitter(alpha=0.1)
    km.fit(a["Death"], a["month"])
    b = km.survival_function_.values.reshape(1,-1)[0]
    ax.line(x=km.survival_function_.index,y=b, label=label, c=c, ls=ls, **kparam)
    if c!='#999999':
        for i in [0,10,20,30,]:
            asum=(a['Death']>i).sum()
            ax.text(i,.2+offset,str(asum),ha='center',color=c)


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
    "Sex (male)",
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
    figwidth='8.4cm',
    share=False,
    gridalpha=.5,
    gridcolor='#dddddd',
    gridstyle='-', #(0,(1,5)),
    gridlinewidth=.3,
    left=None,
    right=None,
    yticks=None,
    yticklen=0,
    labelpad='.25em',
    xloc=('axes',-.025),
    yloc='left',
    refaspect=.73,
    wpad='-.6em',
    xminorlocator='minor',
    yminorlocator='minor',
    xminorlocator_kw={'n': 2},
    yminorlocator_kw={'n': 2},

)
ax = axs[0]
x['z']=1
blues=(x["exp(coef) lower 95%"]<1)
reds=(x["exp(coef) upper 95%"]>1)
lparam=dict(snap=True,lw=.75,alpha=.9,capstyle='round')
yticks=np.array(yticks)
ax.hlines(
    x1=x.loc[blues,"exp(coef) lower 95%"].values,
    x2=x.loc[blues,'z'].values,
    y=yticks[blues.values],
    color="#121212",
    **lparam,

)
ax.hlines(
    x1=x.loc[reds,'z'].values,
    x2=x.loc[reds,"exp(coef) upper 95%"].values,
    y=yticks[reds.values],
    color="#121212",
    **lparam,
)
ax.axvline(1, color="#121212", snap=True, alpha=0.6, lw=.6,solid_capstyle='round')
ax.scatterx(
    x=x["exp(coef)"],
    y=yticks,
    c="#909090",
    s=12,
    marker="D",
    snap=True,
    alpha=0.9,
    mew=0,
)
ylim = list(ax.get_ylim())
ylim[1] += .0


ax.format(
    xlim=(-0, 4.5),
    grid=False,
    ylim=ylim,
    xlabel="Hazard ratio",
    #title="Univariate Cox analysis",
    titleloc="left",
    yloc="none",
    xticks=[0,1,2,3,4],
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


panel = ax.panel("left", share=True, pad='-1.1em')

for i, v in enumerate(x.index):
    sigs = ""
    if x.loc[v, "p"] < 0.05:
        sigs = sigstar(x.loc[v, "p"].astype(float))
    tmp = vdict_r[v]
    if meta[tmp].max() == 1:
        panel.text(
            0.16,
            yticks[i],
            s="(n={}) ".format(meta[tmp].sum().astype(int)) + str(v) + sigs,
            va="center",
            ha="left",
        )  # ,fontsize=7)
    else:
        panel.text(
            0.16, yticks[i], s=str(v) + sigs, va="center", ha="left"
        )  # ,fontsize=7)

panel.format(
    xlim=(-0.02, 2),
    ylim=ax.get_ylim(),
    grid=False,
    xloc="none",
    yloc="none",
)
tparam = dict(
    #weight="demi",
    fontstyle='italic',
    ha="left",
    va="center",
)
panel.text(0, 20.375, s="Microbiome", **tparam)
panel.text(0, 4.375, s="Antibacterials", **tparam)
panel.text(0, 12.375, s="Comorbidities", **tparam)
panel.text(0, 16.375, s="Demographics", **tparam)

ax = axs[1]

ax.format(
    yloc=("axes", -0.07),
    ylim=(0.2, 1.005),
    xlim=(0, 35),
    yticks=np.array([.2,.4,.6,.8,1]),
    ticklen=2,
    yformatter="{x:.0%}",
    xlabel="Days from admission",
    grid=True,
    #title="Kaplan-Meier curves",
    ylabel="Survival",
    ylabelpad=.2,
    yticklabelpad=.2,
    xticks=[0,10,20,30,],

)
for c,v,l,o in zip([blue,red],[(meta['M']==0),(meta['M']==1)],['MMI < 0','MMI > 0'],[.07,.01]):
    plotkm2(ax,v,c=c,label=l,offset=o)
#plotkm2(ax,meta.index,) #ls=(0,(4.5,.85)),)
leg = ax.legend(
    loc='upper right',
    bbox_to_anchor=(1.02,.96),
    ncols=1,
    borderpad=.5,
    borderaxespad=.0,
    handlelength=1.15,
    fontsize=5,
    #frameon=True,
    framealpha=1,
    edgecolor=(1,1,1,0),
)
# drawmedian(ax,)

# %%
fig.savefig("../plots/forest.pdf", dpi=1200, transparent=True, )
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
