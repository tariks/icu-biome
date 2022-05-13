# %%
%matplotlib inline
import proplot as pplt
import numpy as np
import pandas as pd

meta = pd.read_csv("../meta/meta52_current.csv", index_col=0)
x = pd.read_csv('../survival/univariate.csv',index_col=0)
vdict=pd.read_csv('../meta/vdict.csv',index_col=0,header=None)
vdict
# %%
x=x.loc[x.index[::-1]]
x

# %%
fig,axs=pplt.subplots(nrows=1,ncols=3)
ax=axs[1]
ax.hlines(x1=x['exp(coef) lower 95%'], x2=x['exp(coef) upper 95%'],color='#121212')
ax.axvline(1,color='#121212')

# %%
