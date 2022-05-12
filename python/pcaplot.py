# %%
%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
import cmasher as cmr
from deicode.rpca import auto_rpca as rpca
from biom import Table

plt.style.use("../lance.txt")
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


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    a = np.cos(phi) * rho
    b = np.sin(phi) * rho
    ha = "right" if a < 0 else "left"
    va = "top" if b < 0 else "bottom"
    if abs(a) > abs(b):
        va = "center"
    else:
        ha = "center"
    return a, b, ha, va


# %%
meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
table = pd.read_csv("../feature_tables/52N_genus_t80.csv", index_col=0)
T = Table(table.T.values, table.columns, table.index)
mbal = meta["mmi"]
# mbal[mbal>0]=1.5
# mbal[mbal<=0]=-1.5
norm = plt.cm.colors.CenteredNorm(vcenter=0, halfrange=2, clip=True)
norm = plt.cm.colors.Normalize(vmin=-3, vmax=3, clip=True)

arrow = dict(
    width=0.5,
    headwidth=3,
    headlength=3,
    color="#810004",
    alpha=0.5,
)
cmap = plt.get_cmap("coolwarm")
cmap = cmr.apple_r
# %%
oord, dis = rpca(T, min_feature_frequency=30, min_feature_count=100)
pc, L, ex = oord.samples, oord.features, oord.proportion_explained
top = L["PC1"] ** 2 + L["PC2"] ** 2
top = top.sort_values()[::-1]
print(top[:10], top.shape[0])
a = L["PC1"][L["PC1"].abs().sort_values()[::-1].index]
b = L["PC2"][L["PC2"].abs().sort_values()[::-1].index]
c = cmap(norm(mbal[pc.index]))
# c=plt.cm.get_cmap('coolwarm')(norm(mbal[pc.index]))
topa, topb = a.index[0], b.index[0]
fig, ax = plt.subplots()
ax.scatter(
    pc["PC1"],
    pc["PC2"],
    c=c,
    marker="o",
    s=20,
    edgecolors="#bbbbbb",
    lw=0.5,
    snap=True,
    aa=True,
    alpha=0.7,
)
# for tax in [topa,topb]:
for tax in top.index[:5]:
    x, y = 0.52 * a[tax], 0.52 * b[tax]
    ax.annotate(
        "",
        (x, y),
        xytext=(0, 0),
        arrowprops=arrow,
    )
    rho, phi = cart2pol(x, y)
    x, y, ha, va = pol2cart(rho + 0.02, phi)
    ax.text(x, y, tax, ha=ha, va=va, fontsize=5)

ax.set_xlabel("PC1, expl. var. ={:.2%}".format(ex[0]), loc="right")
ax.set_ylabel("PC2, expl. var. ={:.2%}".format(ex[1]), loc="top")
ax.set_title("Robust Aitchison PCA", loc="left")
ax.margins(0.07)
limx, limy = ax.get_xlim(), ax.get_ylim()
ax.set_ylim(min(limy[0], 0.55 * b[topb] - 0.1), max(limy[1], 0.55 * b[topb] + 0.025))
ax.set_xlim(min(limx[0], 0.55 * a[topa] - 0.2), max(limx[1], 0.55 * a[topa] + 0.2))
cbar = plt.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap),
    ax=ax,
    fraction=0.03,
    location="right",
    aspect=25,
    shrink=0.6,
    pad=0.01,
)
cbar.outline.set(lw=0)
plt.savefig("../plots/pca_genus80.pdf", dpi=300, transparent=True, bbox_inches="tight")
# %%
fig, axs = plt.subplots(2, 1, sharex=False)
fig.set(figheight=6, figwidth=2)
fig.subplots_adjust(hspace=0.2)
a = a[:20].sort_values()[::-1]
b = b[:20].sort_values()[::-1]
ax.margins(0.05)
ax = axs[0]
ax.set_xlabel("Magnitude of projection on PC1")
ax.set_title("Top 20 taxa by PC loading")
ax.barh(y=a.index, width=a)
ax = axs[1]
ax.barh(y=b.index, width=b)
ax.set_xlabel("Magnitude of projection on PC2")
ax.margins(0.05)
plt.savefig(
    "../plots/pca_loadings_genus80.pdf", dpi=300, transparent=True, bbox_inches="tight"
)
# %%
meta["M"] = 0
meta.loc[meta["mbal"] > 1.9, "M"] = 1
meta.loc[meta["mbal"] > 4, "M"] = 2
meta.groupby("M").describe()["Death"]
# %%
fig, ax = plt.subplots()
param = dict(s=3, alpha=0.9, edgecolor="#121212", linewidth=0.25, marker="d")
c = sns.color_palette(["#333333", (0.867, 0.667, 0.200, 1.000)])
g = sns.swarmplot(
    data=meta,
    x="M",
    y="Death",
    palette=c,
    hue="month",
    hue_order=[0, 1],
    ax=ax,
    **param
)

fig.set(figheight=1.5, figwidth=1.8)

ax.set_xlabel("Microbial mortality index")
ax.set_xticklabels(["< 0", "< 1.9", "> 4"])
ax.set_ylim(-10, 380)
ax.set_yticks([0, 90, 180, 270, 360])
ax.set_ylabel("Death-free days")
ax.get_legend().remove()
h, l = ax.get_legend_handles_labels()
h = [plt.scatter(1, 1, fc="#333333", **param)]
h.append(plt.scatter(1, 1, fc=(0.867, 0.667, 0.200, 1.000), **param))
ax.legend(
    h,
    ["0", "1"],
    loc="upper right",
    bbox_to_anchor=[1, 1.025],
    title="28-day mortality",
    handletextpad=0.5,
    borderaxespad=0,
    borderpad=0,
    markerscale=2,
    title_fontsize=6,
)
# ax.set_xlim(-.5,2.5)
plt.savefig("../plots/swarmplot.pdf", dpi=300, transparent=True, bbox_inches="tight")
