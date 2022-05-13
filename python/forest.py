# %%
%matplotlib inline
import proplot as pplt
import numpy as np
import pandas as pd

meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
x = pd.read_csv("../survival/univariate.csv", index_col=0)
vdict = pd.read_csv("../meta/vdict.csv", index_col=0, header=None)
vdict
# %%
x = x.loc[x.index[::-1]]
x.index=[vdict.loc[i].values[0] for i in x.index]
idx = ['Ciprofloxacin','Ceftriaxone','Piperacillin-tazobactam','Vancomycin']
idx += ['Diabetes','Sepsis-3','Immunosuppressed','Malignancy','ARDS','APACHE II score']
idx += ['Age','Gender (male)']
idx += ['MMI > 0','MMI']
x=x.loc[idx]
x
# %%
fig, axs = pplt.subplots(nrows=1, ncols=2,journal='ams3')
fig.format(fontsize=6,
           fontfamily='sans',
           linewidth=.6,
           yticks=pplt.Locator('null'),
           )
ax = axs[0,0]
ax.hlines(x1=x["exp(coef) lower 95%"], x2=x["exp(coef) upper 95%"],
                 color="#121212",snap=True,aa=True,lw=1)
ax.axvline(1, color="#121212",snap=True,aa=True,lw=.6)
ax.scatterx(x=x["exp(coef)"],c='#121212',s=x['model LL'],marker='s',smin=10,smax=20,snap=True,aa=True)
panel=ax.panel('left')
for i,v in enumerate(x.index):
    panel.text(0,i-1,s=str(v))
panel.format(
           yticklabels=pplt.Formatter('null'),)
print(fig.get_figwidth(),fig.get_figheight())


# %%
'''
nat1 = 3.5 1.55
nat2 = 7.2 2.79
agu2 = 7.48 4.52
pnas2 = 4.49 1.88
aaas2 = 4.72 1.96
ams3 = 5.5 2.2
ams2 = 4.5 1.89
'''

