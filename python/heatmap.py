# %%
%matplotlib inline
import numpy as np
import pandas as pd
import proplot as pplt
from skbio.stats.composition import clr,multiplicative_replacement
from scipy.cluster.hierarchy import leaves_list, linkage

pplt.rc.load('./.proplotrc')
pplt.config_inline_backend('retina')
pplt.rc.fontfamily = "Source Sans Pro"
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
    yloc="none",
    xrotation=90,
    linewidth=0,
    yticklabels=[],
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
    "lw": 0,
    "ticklabelsize": 7,
    "ticklabels": "{x: .1f}",
    "space": 0.2,
    "tickminor": False,
}

link = linkage(
    corr,
    optimal_ordering=True,
    #metric="canberra",
    #metric="braycurtis",
    #method="complete",
    #method="average",
    #method="single",
    method='ward',
)
order = leaves_list(link)
heat = ax.heatmap(
    corr.iloc[order][cols],
    colorbar="l",
    colorbar_kw=cbk,
    cmap="Div",
    aspect="auto",
    snap=True,
)
ax.hlines(y=np.arange(-.5,31,1),x1=-.5,x2=4.5,lw=.5,ec='#ffffff',snap=True,alpha=.75)
ax.vlines(x=[-.5,1.5,2.5,4.5],y1=-.5,y2=30.5,lw=4.8,color='#ffffff',snap=True)
ax.vlines(x=[1.525,2.475],y1=-.5,y2=30.5,lw=.3,color='#000000',snap=True)
labs = [i.replace('[','') for i in corr.index[order]]
ax.format(yloc='right',yticklabels=labs)
newlabs=[]
for j,i in enumerate(labs):
    c='#121212'
    if i in aerobes:
        c='blueviolet'
    if i in gpac:
        c='orangered'
    ax.text(4.5,j,i,c=c,va='center',fontsize=6)
    newlabs.append(i)
ax.set_yticklabels([])
ax.text(.5,30.75,'Numerator',ha='center',va='bottom',fontsize=6.5)
ax.text(3.5,30.75,'Denominator',ha='center',va='bottom',fontsize=6.5)


fig.set(figheight=3.4,figwidth=2.85)
#ax.invert_yaxis()
# %%
fig.savefig(
    "../plots/taxaheatmap.pdf",
    dpi=600,
    transparent=True,
    bbox_inches='tight'
)


# %%
corr.to_csv('../plots/heatmap.csv')
# %%
import pandas as pd

# %%
meta=pd.read_csv('../meta/meta52_current.csv',index_col=0)
meta.columns
# %%
cols = [
       'Age',
       'Gender',
       'Race',
       'APACHE',
       'PF',
       'ARDS',
       'Sepsis3',
       'Shock', 
       'malignancy',
       'diabetes', 
       'immunosuppressed',
       'COPD',
       'Zithromax',
       'Cefepime', 
       'Ceftriaxone', 
       'Cipro',
       'Zosyn', 
       'Vancomycin',
       'mmi',
       'M', 
       'month', 
       'year',
       'Death',
]
m = meta[cols]
m.loc[m['Race']==0,'Race'] = 5
m.loc[m['Race']<5,'Race'] = 0
m.loc[m['Race']>0,'Race'] = 1
# %%
nums = ['Age','PF','Death','APACHE','mmi']
cats = [i for i in cols if i not in nums]
vdict = {i: i for i in cols}
vdict["mmi"] = "MMI"
vdict["M"] = "MMI > 0"
vdict["month"] = "Death < 28 days"
vdict["year"] = "Death < 1 year"
vdict["Cipro"] = "Ciprofloxacin"
vdict["Zithromax"] = "Azithromycin"
vdict["Zosyn"] = "Piperacillin-tazobactam"
vdict["immunosuppressed"] = "Immunosuppressed"
vdict["diabetes"] = "Diabetes"
vdict["malignancy"] = "Malignancy"
vdict["Death"] = "Death-free days"
vdict["APACHE"] = "APACHE II score"
vdict["Sepsis3"] = "Sepsis-3"
vdict["Gender"] = "Gender (male)"
vdict["Shock"] = "Septic shock"
vdict["PF"] = "PF ratio"
vdict["Race"] = "Race, caucasian"
# %%
x = pd.Series()
x['N'] = 52
for i in cols:
    name = vdict.get(i)
    if i in nums:
        median = m[i].median()
        iqr = m[i].quantile(.75) - m[i].quantile(.25)
        x[name] = '{:.1f} [{:.1f}]'.format(median,iqr)
    else:
        pct = m[i].sum()/52
        s = m[i].sum().astype(int)
        x[name] = '{} ({:.0%})'.format(s,pct)
x

# %%
